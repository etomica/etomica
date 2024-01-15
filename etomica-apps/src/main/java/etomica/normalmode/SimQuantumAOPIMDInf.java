/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.normalmode;

import etomica.action.activity.ActivityIntegrate;
import etomica.atom.AtomType;
import etomica.atom.DiameterHashByType;
import etomica.atom.IAtom;
import etomica.box.Box;
import etomica.config.ConformationLinear;
import etomica.data.*;
import etomica.data.types.DataGroup;
import etomica.graphics.ColorScheme;
import etomica.graphics.DisplayBox;
import etomica.graphics.DisplayTextBox;
import etomica.graphics.SimulationGraphic;
import etomica.integrator.Integrator;
import etomica.integrator.IntegratorLangevin;
import etomica.integrator.IntegratorMD;
import etomica.potential.*;
import etomica.potential.compute.PotentialCompute;
import etomica.potential.compute.PotentialComputeAggregate;
import etomica.potential.compute.PotentialComputeField;
import etomica.simulation.Simulation;
import etomica.space.BoundaryRectangularNonperiodic;
import etomica.space.Space;
import etomica.space1d.Space1D;
import etomica.species.SpeciesBuilder;
import etomica.species.SpeciesGeneral;
import etomica.util.ParameterBase;
import etomica.util.ParseArgs;

import java.awt.*;
import java.util.ArrayList;
import java.util.List;

public class SimQuantumAOPIMDInf extends Simulation {

    public PotentialComputeField pcP1harm, pcP1ah, pcP1EnTIA;
    public final IntegratorMD integrator;
    public final PotentialMasterBonding pmBonding;
    public final Box box;
    public final PotentialComputeAggregate pmAgg;
    public P1AnharmonicTIA p1ahUeff;
    public int dim, nBeads;
    public MCMoveHOReal2 moveStageSimple, moveStageEC;

    public SimQuantumAOPIMDInf(Space space, MoveChoice coordType, double mass, double timeStep, double gammaLangevin, int nBeads, double temperature, double k4, double omega, double mOmegaF2, double mOmegaH2, boolean isTIA, double hbar) {
        super(space);

        SpeciesGeneral species = new SpeciesBuilder(space)
                .setDynamic(true)
                .addCount(AtomType.simple("A", mass / nBeads), nBeads)
                .withConformation(new ConformationLinear(space, 0))
                .build();
        addSpecies(species);
        this.nBeads = nBeads;
        box = this.makeBox(new BoundaryRectangularNonperiodic(space));
        box.setNMolecules(species, 1);
        //pm2 that uses the full PI potential, for data collection
        //spring P2 part (x_i-x_{i+1})^2
        pmBonding = new PotentialMasterBonding(getSpeciesManager(), box);
        this.dim = space.D();

        P2Harmonic p2Bond = new P2Harmonic(mOmegaF2, 0);
        List<int[]> pairs = new ArrayList<>();
        for (int i=0; i<nBeads; i++) {
            int[] p = new int[]{i,(i+1)%nBeads};
            pairs.add(p);
        }
        pmBonding.setBondingPotentialPair(species, p2Bond, pairs);
        pcP1harm = new PotentialComputeField(getSpeciesManager(), box);
        pcP1ah = new PotentialComputeField(getSpeciesManager(), box);

        IPotential1 p1harm, p1ah;

        if (isTIA){
//            double facUeff = 1.0;
//            p1ahUeff = new P1AnharmonicTIA(space, 1, k4, nBeads, mass*omegaN2, facUeff);
//            pcP1.setFieldPotential(species.getLeafType(), p1ahUeff);
        } else {
            p1harm = new P1Anharmonic(space, mOmegaH2, 0);
            pcP1harm.setFieldPotential(species.getLeafType(), p1harm);
            p1ah = new P1Anharmonic(space, 0, k4/nBeads);
            pcP1ah.setFieldPotential(species.getLeafType(), p1ah);
        }

        PotentialComputeAggregate.localStorageDefault = true;
        pmAgg = new PotentialComputeAggregate(pmBonding, pcP1harm, pcP1ah);

//        double facEn = 3.0;
//        P1AnharmonicTIA p1ahEn = new P1AnharmonicTIA(space, 1, k4, nBeads, mass*omegaN2, facEn);
//        pcP1EnTIA = new PotentialComputeField(getSpeciesManager(), box);
//        pcP1EnTIA.setFieldPotential(species.getLeafType(), p1ahEn);

        double omegaSample = omega;

        if (coordType == MoveChoice.Real) {
            integrator = new IntegratorLangevin(pmAgg, random, timeStep, temperature, box, gammaLangevin);
        } else if (coordType == MoveChoice.NM) {
            integrator = new IntegratorLangevinPINMInf(pmAgg, random, timeStep, temperature, box, gammaLangevin, hbar, 0, omegaSample);
        } else if (coordType == MoveChoice.NMEC) {
            integrator = new IntegratorLangevinPINMInf(pmAgg, random, timeStep, temperature, box, gammaLangevin, hbar, omega, omegaSample);
        } else if (coordType == MoveChoice.Stage) {
            integrator = new IntegratorLangevinPIInf(pmAgg, random, timeStep, temperature, box, gammaLangevin, hbar, omega, 0);
        } else { //StageEC -- default
            integrator = new IntegratorLangevinPIInf(pmAgg, random, timeStep, temperature, box, gammaLangevin,  hbar, omega, omegaSample);
        }

        moveStageEC = new MCMoveHOReal2(space, pmAgg, random, temperature, omega*omega, box, hbar);
        moveStageSimple = new MCMoveHOReal2(space, pmAgg, random, temperature, 0, box, hbar);
        integrator.setThermostatNoDrift(false);
        integrator.setIsothermal(true);
    }

    public Integrator getIntegrator() {
        return integrator;
    }

    public static void main(String[] args) {
        final long startTime = System.currentTimeMillis();
        OctaneParams params = new OctaneParams();
        if (args.length > 0) {
            ParseArgs.doParseArgs(params, args);
        } else {
            // custom parameters
            params.hbar = 1;
            params.steps = 100000;
            params.temperature = 1;
            params.omega = 1;
            params.k4 = 0;
//            params.coordType = MoveChoice.Real;
//            params.coordType = MoveChoice.NMEC;
            params.coordType = MoveChoice.StageEC;
        }

        int nShifts = params.nShifts;
        double mass = params.mass;
        int nBeads = params.nBeads;
        double temperature = params.temperature;
        double beta = 1/temperature;
        double hbar = params.hbar;
        double omega = params.omega;
        double k4 = params.k4;
        boolean isGraphic = params.isGraphic;
        long steps = params.steps;
        long stepsEq = steps/10;
        boolean isTIA = params.isTIA;
        boolean zerok0 = params.zerok0;
        boolean onlyCentroid = params.onlyCentroid;
        MoveChoice coordType = params.coordType;
        double omega2 = omega*omega;

        double x = hbar*omega/temperature;
        if (nBeads == -1){
            nBeads = (int) (20*x); //20*x and 30*x are good for HO and AO, resp.
        }
        double alpha = beta*hbar*omega/nBeads;
        double mOmegaF2 = nBeads == 1 ? 0 : mass*omega2/nBeads/alpha/Math.sinh(alpha);
        double mOmegaH2 = 2*mass*omega2*Math.tanh(alpha/2)/nBeads/alpha;
        double omegaRing = Math.sqrt(mOmegaH2*nBeads/mass);
        double omegaBead = Math.sqrt(mOmegaF2*nBeads/mass);

        double omegaN = Math.sqrt(nBeads)*temperature/hbar;
        double timeStep = params.timeStep;
        if (timeStep == -1) {
            double c = 0.1;
            if (params.coordType == MoveChoice.Real) {
                timeStep = c/omegaBead;// mi=m/n in real space, so m wn^2 = (m/n)*(n wn^2)==> dt~1/[wn sqrt(n)]
            } else if (params.coordType == MoveChoice.NM) {
                double s = omega2/nBeads;
                timeStep = c*Math.sqrt(s)/omega; // which is 1/sqrt(n)
            } else if (params.coordType == MoveChoice.Stage) {
                double s = omega2/omegaN/omegaN*(1 + 1.0/12.0*(nBeads*nBeads-1.0)/nBeads);
                timeStep = c*Math.sqrt(s)/omega; // for large n, timeStep ~ hbar/T
            } else {
                timeStep = c/omegaRing;
            }
        }

//        if (isTIA){
//            omega2 = omega2*(1.0 + omega2/12.0/(nBeads*omegaN2));
//        }
//        if (zerok0) omega2 = 0;

        double gammaLangevin = 2*omegaRing;
        final SimQuantumAOPIMDInf sim = new SimQuantumAOPIMDInf(Space1D.getInstance(), coordType, mass, timeStep, gammaLangevin, nBeads, temperature, k4, omega, mOmegaF2, mOmegaH2, isTIA, hbar);
        sim.integrator.reset();
        System.out.println(" PIMDinf-" + coordType);
        System.out.println(" mass: " + mass);
        System.out.println(" T: " + temperature);
        System.out.println(" hbar: " + hbar);
        System.out.println(" w: " + omega);
        System.out.println(" mOmegaF2: " + mOmegaF2  + " , mOmegaH2: " + mOmegaH2);
        System.out.println(" x = beta*hbar*w = " + hbar*omega/temperature);
        System.out.println(" nBeads: " + nBeads);
        System.out.println(" nShifts: "+ nShifts);
        System.out.println(" steps: " +  steps + " stepsEq: " + stepsEq);
        System.out.println(" timestep: " + timeStep);
        System.out.println(" k4: " + k4);
        System.out.println(" isTIA: " + isTIA);
        System.out.println(" gammaLangevin: " + gammaLangevin);

        System.out.println("\n Quantum Harmonic Oscillator Theory");

        DataSourceScalar meterKE = sim.integrator.getMeterKineticEnergy();

//        MeterMSDHO meterMSDHO = new MeterMSDHO(sim.box);
        MeterPIPrimInf meterPrimInf = null;
        MeterPIVirInf meterVirInf = null;
        MeterPICentVirInf meterCentVirInf = null;
        MeterPIVirHarmStagingInf meterVirHarmStagingInf = null;
        MeterPIHMAcInf meterHMAcInf = null;
        MeterPIHMAInf meterNMECInf = null;
        MeterPIHMAInf meterNMSimpleInf = null;
        MeterPIHMAReal2Inf meterStageECInf = null;
        MeterPIHMAReal2Inf meterStageSimpleInf = null;
        if (isTIA){
//            meterPrim = new MeterPIPrim(sim.pmBonding, sim.pcP1EnTIA, nBeads, sim.betaN);
//            meterVir = new MeterPIVirTIA(sim.pcP1EnTIA, sim.pcP1, sim.betaN, nBeads, sim.box);
//            meterCentVir = new MeterPICentVirTIA(sim.pcP1EnTIA, sim.pcP1, sim.betaN, nBeads, sim.box);
            meterHMAcInf = null;
//            meterHMA = new MeterPIHMATIA(sim.pmBonding, sim.pcP1EnTIA, sim.pcP1, sim.betaN, nBeads, omega2, sim.box);
//            meterHMAcent = null;
        } else {
            meterPrimInf = new MeterPIPrimInf(sim.pmBonding, sim.pcP1harm, sim.pcP1ah, nBeads, temperature, sim.box, omega, hbar);
            meterVirInf = new MeterPIVirInf(sim.pcP1harm, sim.pcP1ah, temperature, sim.box,  nBeads, omega, hbar);
            meterCentVirInf = new MeterPICentVirInf(sim.pcP1harm, sim.pcP1ah, temperature, sim.box, nBeads, omega, hbar);
            meterVirHarmStagingInf = new MeterPIVirHarmStagingInf(sim.pcP1ah, temperature, sim.box,  nBeads, omega, hbar);
            meterHMAcInf = new MeterPIHMAcInf(sim.pcP1harm, sim.pcP1ah, temperature, sim.box, nBeads, omega, hbar);
            meterNMSimpleInf = new MeterPIHMAInf(sim.pmBonding, sim.pcP1harm, sim.pcP1ah, temperature, nBeads, omega, 0, sim.box, hbar);
            meterNMECInf = new MeterPIHMAInf(sim.pmBonding, sim.pcP1harm, sim.pcP1ah, temperature, nBeads, omega, omega, sim.box, hbar);
            meterStageSimpleInf = new MeterPIHMAReal2Inf(sim.pmBonding, sim.pcP1harm, sim.pcP1ah, temperature, nBeads, omega, 0, sim.box, hbar);
            meterStageSimpleInf.setNumShifts(nShifts);
            meterStageECInf = new MeterPIHMAReal2Inf(sim.pmBonding, sim.pcP1harm, sim.pcP1ah, temperature, nBeads, omega, omega, sim.box, hbar);
            meterStageECInf.setNumShifts(nShifts);
        }

        if (isGraphic) {
            sim.getController().addActivity(new ActivityIntegrate(sim.integrator));
            SimulationGraphic simGraphic = new SimulationGraphic(sim, SimulationGraphic.TABBED_PANE);
            simGraphic.setPaintInterval(sim.box, 1);
            int finalNBeads = nBeads;
            ColorScheme colorScheme = new ColorScheme() {
                protected Color[] allColors;

                public Color getAtomColor(IAtom a) {
                    if (allColors == null) {
                        allColors = new Color[768];
                        for (int i = 0; i < 256; i++) {
                            allColors[i] = new Color(255 - i, i, 0);
                        }
                        for (int i = 0; i < 256; i++) {
                            allColors[i + 256] = new Color(0, 255 - i, i);
                        }
                        for (int i = 0; i < 256; i++) {
                            allColors[i + 512] = new Color(i, 0, 255 - i);
                        }
                    }
                    return allColors[(768 * a.getIndex() / (finalNBeads))];
                }
            };

            simGraphic.getDisplayBox(sim.box).setColorScheme(colorScheme);
            ((DiameterHashByType) ((DisplayBox) simGraphic.displayList().getFirst()).getDiameterHash()).setDiameter(sim.species().getAtomType(0), 1.0);
            DisplayTextBox timer = new DisplayTextBox();
            DataSourceCountTime counter = new DataSourceCountTime(sim.integrator);
            DataPumpListener counterPump = new DataPumpListener(counter, timer, 100);
            sim.integrator.getEventManager().addListener(counterPump);
            simGraphic.getPanel().controlPanel.add(timer.graphic());
            simGraphic.makeAndDisplayFrame("PIMD - "+coordType);
            return;
        }

        System.out.flush();
        // equilibration
        sim.getController().runActivityBlocking(new ActivityIntegrate(sim.integrator, stepsEq));
        System.out.println("\n equilibration finished");

        int interval = 10;
        int blocks = 100;
        long blockSize = steps / (interval * blocks);
        if (blockSize == 0) blockSize = 1;
        System.out.println(" numBlocks: " + blocks + " blocksize: " + blockSize + " interval: " + interval);

//        AccumulatorAverageFixed accumulatorMSD = new AccumulatorAverageFixed(blockSize);
//        DataPumpListener accumulatorPumpMSD = new DataPumpListener(meterMSDHO, accumulatorMSD, interval);
//        sim.integrator.getEventManager().addListener(accumulatorPumpMSD);


        AccumulatorAverageCovariance accumulatorKE = new AccumulatorAverageCovariance(blockSize);
        DataPumpListener accumulatorPumpKE = new DataPumpListener(meterKE, accumulatorKE, interval);
        sim.integrator.getEventManager().addListener(accumulatorPumpKE);


        //1 Primitive
        AccumulatorAverageCovariance accumulatorPrimInf = new AccumulatorAverageCovariance(blockSize);
        DataPumpListener accumulatorPumpPrimInf = new DataPumpListener(meterPrimInf, accumulatorPrimInf, interval);
        sim.integrator.getEventManager().addListener(accumulatorPumpPrimInf);

        //2 Virial
        AccumulatorAverageCovariance accumulatorVirInf = new AccumulatorAverageCovariance(blockSize);
        if (meterVirInf != null) {
            DataPumpListener accumulatorPumpVirInf = new DataPumpListener(meterVirInf, accumulatorVirInf, interval);
            sim.integrator.getEventManager().addListener(accumulatorPumpVirInf);
        }

        AccumulatorAverageCovariance accumulatorVirHarmStagingInf = new AccumulatorAverageCovariance(blockSize);
        if (meterVirHarmStagingInf != null) {
            DataPumpListener accumulatorPumpVirHarmStagingInf = new DataPumpListener(meterVirHarmStagingInf, accumulatorVirHarmStagingInf, interval);
            sim.integrator.getEventManager().addListener(accumulatorPumpVirHarmStagingInf);
        }

        //3 Centroid Virial
        AccumulatorAverageCovariance accumulatorCentVirInf = new AccumulatorAverageCovariance(blockSize);
        if (meterCentVirInf != null) {
            DataPumpListener accumulatorPumpCentVirInf = new DataPumpListener(meterCentVirInf, accumulatorCentVirInf, interval);
            sim.integrator.getEventManager().addListener(accumulatorPumpCentVirInf);
        }

         //4 HMAc (CLassical EC)
        AccumulatorAverageCovariance accumulatorHMAcInf = new AccumulatorAverageCovariance(blockSize);
        if (meterHMAcInf != null) {
            DataPumpListener accumulatorPumpHMAcInf = new DataPumpListener(meterHMAcInf, accumulatorHMAcInf, interval);
            sim.integrator.getEventManager().addListener(accumulatorPumpHMAcInf);
        }

        //5 HMAq (Quantum EC)
        AccumulatorAverageCovariance accumulatorNMSimpleInf = new AccumulatorAverageCovariance(blockSize);
        if (meterNMSimpleInf != null) {
            DataPumpListener accumulatorPumpHMAsimpleInf = new DataPumpListener(meterNMSimpleInf, accumulatorNMSimpleInf, interval);
            sim.integrator.getEventManager().addListener(accumulatorPumpHMAsimpleInf);
        }

        AccumulatorAverageCovariance accumulatorNMECInf = new AccumulatorAverageCovariance(blockSize);
        if (meterNMECInf != null) {
            DataPumpListener accumulatorPumpHMAInf = new DataPumpListener(meterNMECInf, accumulatorNMECInf, interval);
            sim.integrator.getEventManager().addListener(accumulatorPumpHMAInf);
        }

        AccumulatorAverageCovariance accumulatorStageSimpleInf = new AccumulatorAverageCovariance(blockSize);
        if (meterStageSimpleInf != null) {
            DataPumpListener pumpStagesimpleInf = new DataPumpListener(meterStageSimpleInf, accumulatorStageSimpleInf, interval);
            sim.integrator.getEventManager().addListener(pumpStagesimpleInf);
        }

        AccumulatorAverageCovariance accumulatorStageECInf = new AccumulatorAverageCovariance(blockSize);
        if (meterStageECInf != null) {
            DataPumpListener pumpStageECInf = new DataPumpListener(meterStageECInf, accumulatorStageECInf, interval);
            sim.integrator.getEventManager().addListener(pumpStageECInf);
        }



        //run
        sim.integrator.resetStepCount();
        sim.getController().runActivityBlocking(new ActivityIntegrate(sim.integrator, steps));

        //MSD
//        DataGroup dataMSD = (DataGroup)accumulatorMSD.getData();
//        IData dataMSDAvg = dataMSD.getData(accumulatorMSD.AVERAGE.index);
//        IData dataMSDErr = dataMSD.getData(accumulatorMSD.ERROR.index);
//        IData dataMSDCorrelation = dataMSD.getData(accumulatorMSD.BLOCK_CORRELATION.index);
//        double avgMSD = dataMSDAvg.getValue(0);
//        double errMSD = dataMSDErr.getValue(0);
//        double corMSD = dataMSDCorrelation.getValue(0);

        //T
        DataGroup dataKE = (DataGroup) accumulatorKE.getData();
        double avgKE = dataKE.getValue(accumulatorKE.AVERAGE.index);
        double errKE = dataKE.getValue(accumulatorKE.ERROR.index);
        double corKE = dataKE.getValue(accumulatorKE.BLOCK_CORRELATION.index);
        System.out.println("\n T_measured: " + avgKE / (0.5*sim.dim* nBeads) + " +/- " + errKE / (0.5*sim.dim*nBeads) + " cor: " + corKE);

        double kB_beta2 = beta*beta;
        double varX0, varX1, corX0X1;

        //1 Prim
        DataGroup dataPrimInf = (DataGroup) accumulatorPrimInf.getData();
        IData dataAvgPrimInf = dataPrimInf.getData(accumulatorPrimInf.AVERAGE.index);
        IData dataErrPrimInf = dataPrimInf.getData(accumulatorPrimInf.ERROR.index);
        IData dataCorPrimInf = dataPrimInf.getData(accumulatorPrimInf.BLOCK_CORRELATION.index);
        IData dataCovPrimInf = dataPrimInf.getData(accumulatorPrimInf.COVARIANCE.index);
        double avgEnPrimInf = dataAvgPrimInf.getValue(0);
        double errEnPrimInf = dataErrPrimInf.getValue(0);
        double corEnPrimInf = dataCorPrimInf.getValue(0);
        System.out.println("\n En_prim:          " + avgEnPrimInf + "   err: " + errEnPrimInf + " cor: " + corEnPrimInf);
        double CvnPrimInf = kB_beta2*(dataAvgPrimInf.getValue(1) - avgEnPrimInf*avgEnPrimInf);
        varX0 = errEnPrimInf*errEnPrimInf;
        varX1 = dataErrPrimInf.getValue(1)*dataErrPrimInf.getValue(1);
        corX0X1 = dataCovPrimInf.getValue(1)/Math.sqrt(dataCovPrimInf.getValue(0))/Math.sqrt(dataCovPrimInf.getValue(3));
        double errCvnPrimInf = kB_beta2*Math.sqrt(varX1 + 4.0*avgEnPrimInf*avgEnPrimInf*varX0 - 4*avgEnPrimInf*dataErrPrimInf.getValue(0)*dataErrPrimInf.getValue(1)*corX0X1);


        //2 Vir
        DataGroup dataVirInf = (DataGroup) accumulatorVirInf.getData();
        IData dataAvgVirInf = dataVirInf.getData(accumulatorVirInf.AVERAGE.index);
        IData dataErrVirInf = dataVirInf.getData(accumulatorVirInf.ERROR.index);
        IData dataCorVirInf = dataVirInf.getData(accumulatorVirInf.BLOCK_CORRELATION.index);
        IData dataCovVirInf = dataVirInf.getData(accumulatorVirInf.COVARIANCE.index);
        double avgEnVirInf = dataAvgVirInf.getValue(0);
        double errEnVirInf = dataErrVirInf.getValue(0);
        double corEnVirInf = dataCorVirInf.getValue(0);
        System.out.println(" En_vir:           " + avgEnVirInf + "   err: " + errEnVirInf + " cor: " + corEnVirInf);
        double CvnVirInf = kB_beta2*(dataAvgVirInf.getValue(1) - avgEnVirInf*avgEnVirInf);
        varX0 = errEnVirInf*errEnVirInf;
        varX1 = dataErrVirInf.getValue(1)*dataErrVirInf.getValue(1);
        corX0X1 = dataCovVirInf.getValue(1)/Math.sqrt(dataCovVirInf.getValue(0))/Math.sqrt(dataCovVirInf.getValue(3));
        double errCvnVirInf = kB_beta2*Math.sqrt(varX1 + 4.0*avgEnVirInf*avgEnVirInf*varX0 - 4*avgEnVirInf*dataErrVirInf.getValue(0)*dataErrVirInf.getValue(1)*corX0X1);

        //3 Cent-Vir
        DataGroup dataCentVirInf = (DataGroup) accumulatorCentVirInf.getData();
        IData dataAvgCentVirInf = dataCentVirInf.getData(accumulatorCentVirInf.AVERAGE.index);
        IData dataErrCentVirInf = dataCentVirInf.getData(accumulatorCentVirInf.ERROR.index);
        IData dataCorCentVirInf = dataCentVirInf.getData(accumulatorCentVirInf.BLOCK_CORRELATION.index);
        IData dataCovCentVirInf = dataCentVirInf.getData(accumulatorCentVirInf.COVARIANCE.index);
        double avgEnCentVirInf = dataAvgCentVirInf.getValue(0);
        double errEnCentVirInf = dataErrCentVirInf.getValue(0);
        double corEnCentVirInf = dataCorCentVirInf.getValue(0);
        System.out.println(" En_cvir:          " + avgEnCentVirInf + "   err: " + errEnCentVirInf + " cor: " + corEnCentVirInf);
        double CvnCentVirInf = kB_beta2*(dataAvgCentVirInf.getValue(1) - avgEnCentVirInf*avgEnCentVirInf);
        varX0 = errEnCentVirInf*errEnCentVirInf;
        varX1 = dataErrCentVirInf.getValue(1)*dataErrCentVirInf.getValue(1);
        corX0X1 = dataCovCentVirInf.getValue(1)/Math.sqrt(dataCovCentVirInf.getValue(0))/Math.sqrt(dataCovCentVirInf.getValue(3));
        double errCvnCentVirInf = kB_beta2*Math.sqrt(varX1 + 4.0*avgEnCentVirInf*avgEnCentVirInf*varX0 - 4*avgEnCentVirInf*dataErrCentVirInf.getValue(0)*dataErrCentVirInf.getValue(1)*corX0X1);

        //2' HO Staging Vir
        DataGroup dataVirHarmStagingInf = (DataGroup) accumulatorVirHarmStagingInf.getData();
        IData dataAvgVirHarmStagingInf = dataVirHarmStagingInf.getData(accumulatorVirHarmStagingInf.AVERAGE.index);
        IData dataErrVirHarmStagingInf = dataVirHarmStagingInf.getData(accumulatorVirHarmStagingInf.ERROR.index);
        IData dataCorVirHarmStagingInf = dataVirHarmStagingInf.getData(accumulatorVirHarmStagingInf.BLOCK_CORRELATION.index);
        IData dataCovVirHarmStagingInf = dataVirHarmStagingInf.getData(accumulatorVirHarmStagingInf.COVARIANCE.index);
        double avgEnVirHarmStagingInf = dataAvgVirHarmStagingInf.getValue(0);
        double errEnVirHarmStagingInf = dataErrVirHarmStagingInf.getValue(0);
        double corEnVirHarmStagingInf = dataCorVirHarmStagingInf.getValue(0);
        System.out.println(" En_ho_stage_vir:  " + avgEnVirHarmStagingInf + "   err: " + errEnVirHarmStagingInf + " cor: " + corEnVirHarmStagingInf);

        //4 HMAc
        DataGroup dataHMAcInf = (DataGroup) accumulatorHMAcInf.getData();
        IData dataAvgHMAcInf = dataHMAcInf.getData(accumulatorHMAcInf.AVERAGE.index);
        IData dataErrHMAcInf = dataHMAcInf.getData(accumulatorHMAcInf.ERROR.index);
        IData dataCorHMAcInf = dataHMAcInf.getData(accumulatorHMAcInf.BLOCK_CORRELATION.index);
        IData dataCovHMAcInf = dataHMAcInf.getData(accumulatorHMAcInf.COVARIANCE.index);
        double avgEnHMAcInf = dataAvgHMAcInf.getValue(0) ;
        double errEnHMAcInf = dataErrHMAcInf.getValue(0);
        double corEnHMAcInf = dataCorHMAcInf.getValue(0);
        System.out.println(" En_hmac:          " + avgEnHMAcInf + "   err: " + errEnHMAcInf + " cor: " + corEnHMAcInf);
        double CvnHMAcInf = kB_beta2*(dataAvgHMAcInf.getValue(1) - avgEnHMAcInf*avgEnHMAcInf);
        varX0 = errEnHMAcInf*errEnHMAcInf;
        varX1 = dataErrHMAcInf.getValue(1)*dataErrHMAcInf.getValue(1);
        corX0X1 = dataCovHMAcInf.getValue(1)/Math.sqrt(dataCovHMAcInf.getValue(0))/Math.sqrt(dataCovHMAcInf.getValue(3));
        double errCvnHMAcInf = kB_beta2*Math.sqrt(varX1 + 4.0*avgEnHMAcInf*avgEnHMAcInf*varX0 - 4*avgEnHMAcInf*dataErrHMAcInf.getValue(0)*dataErrHMAcInf.getValue(1)*corX0X1);
        if (errEnHMAcInf < 1e-10){
            errCvnHMAcInf = kB_beta2*Math.sqrt(varX1 + 4.0*avgEnHMAcInf*avgEnHMAcInf*varX0);
        }


        //5 HMA simple NM
        DataGroup dataNMsimpleInf = (DataGroup) accumulatorNMSimpleInf.getData();
        IData dataAvgNMsimpleInf = dataNMsimpleInf.getData(accumulatorNMSimpleInf.AVERAGE.index);
        IData dataErrNMsimpleInf = dataNMsimpleInf.getData(accumulatorNMSimpleInf.ERROR.index);
        IData dataCorNMsimpleInf = dataNMsimpleInf.getData(accumulatorNMSimpleInf.BLOCK_CORRELATION.index);
        IData dataCovNMsimpleInf = dataNMsimpleInf.getData(accumulatorNMSimpleInf.COVARIANCE.index);
        double avgEnNMSimpleInf = dataAvgNMsimpleInf.getValue(0);
        double errEnNMSimpleInf = dataErrNMsimpleInf.getValue(0);
        double corEnNMSimpleInf = dataCorNMsimpleInf.getValue(0);
        System.out.println(" En_nm_simple:     " + avgEnNMSimpleInf + "   err: " + errEnNMSimpleInf + " cor: " + corEnNMSimpleInf);
        double Cvn_nm_simpleInf  = kB_beta2*(dataAvgNMsimpleInf.getValue(1) - avgEnNMSimpleInf*avgEnNMSimpleInf);
        varX0 = errEnNMSimpleInf*errEnNMSimpleInf;
        varX1 = dataErrNMsimpleInf.getValue(1)*dataErrNMsimpleInf.getValue(1);
        corX0X1 = dataCovNMsimpleInf.getValue(1)/Math.sqrt(dataCovNMsimpleInf.getValue(0))/Math.sqrt(dataCovNMsimpleInf.getValue(3));
        double errCvnNMsimpleInf = kB_beta2*Math.sqrt(varX1 + 4.0*avgEnNMSimpleInf*avgEnNMSimpleInf*varX0 - 4*avgEnNMSimpleInf*dataErrNMsimpleInf.getValue(0)*dataErrNMsimpleInf.getValue(1)*corX0X1);

        //6 HMA EC NM
        DataGroup dataNMEC = (DataGroup) accumulatorNMECInf.getData();
        IData dataAvgNMEC = dataNMEC.getData(accumulatorNMECInf.AVERAGE.index);
        IData dataErrNMEC = dataNMEC.getData(accumulatorNMECInf.ERROR.index);
        IData dataCorNMEC = dataNMEC.getData(accumulatorNMECInf.BLOCK_CORRELATION.index);
        IData dataCovNMEC = dataNMEC.getData(accumulatorNMECInf.COVARIANCE.index);
        double avgEnNMEC = dataAvgNMEC.getValue(0);
        double errEnNMEC = dataErrNMEC.getValue(0);
        double corEnNMEC = dataCorNMEC.getValue(0);
        System.out.println(" En_nm_ec:         " + avgEnNMEC + "   err: " + errEnNMEC + " cor: " + corEnNMEC);
        double CvnNMECInf  = kB_beta2*(dataAvgNMEC.getValue(1) - avgEnNMEC*avgEnNMEC);
        varX0 = errEnNMEC*errEnNMEC;
        varX1 = dataErrNMEC.getValue(1)*dataErrNMEC.getValue(1);
        corX0X1 = dataCovNMEC.getValue(1)/Math.sqrt(dataCovNMEC.getValue(0))/Math.sqrt(dataCovNMEC.getValue(3));
        double errCvnNMECInf = kB_beta2*Math.sqrt(varX1 + 4.0*avgEnNMEC*avgEnNMEC*varX0 - 4*avgEnNMEC*errEnNMEC*dataErrNMEC.getValue(1)*corX0X1);
        if (errEnNMEC < 1e-10){
            errCvnNMECInf = kB_beta2*Math.sqrt(varX1 + 4.0*avgEnNMEC*avgEnNMEC*varX0);
        }

        // 7 HMA simple stage
        DataGroup dataStageSimpleInf = (DataGroup) accumulatorStageSimpleInf.getData();
        IData dataAvgStageSimpleInf = dataStageSimpleInf.getData(accumulatorStageSimpleInf.AVERAGE.index);
        IData dataErrStageSimpleInf = dataStageSimpleInf.getData(accumulatorStageSimpleInf.ERROR.index);
        IData dataCorStageSimpleInf = dataStageSimpleInf.getData(accumulatorStageSimpleInf.BLOCK_CORRELATION.index);
        IData dataCovStageSimpleInf = dataStageSimpleInf.getData(accumulatorStageSimpleInf.COVARIANCE.index);
        double avgEnStageSimpleInf = dataAvgStageSimpleInf.getValue(0);
        double errEnStageSimpleInf = dataErrStageSimpleInf.getValue(0);
        double corEnStageSimpleInf = dataCorStageSimpleInf.getValue(0);
        System.out.println(" En_stage_simple:  " + avgEnStageSimpleInf + "   err: " + errEnStageSimpleInf + " cor: " + corEnStageSimpleInf);
        double Cvn_stage_simple  = kB_beta2*(dataAvgStageSimpleInf.getValue(1) - avgEnStageSimpleInf*avgEnStageSimpleInf);
        varX0 = errEnStageSimpleInf*errEnStageSimpleInf;
        varX1 = dataErrStageSimpleInf.getValue(1)*dataErrStageSimpleInf.getValue(1);
        corX0X1 = dataCovStageSimpleInf.getValue(1)/Math.sqrt(dataCovStageSimpleInf.getValue(0))/Math.sqrt(dataCovStageSimpleInf.getValue(3));
        double errCvnStageSimple = kB_beta2*Math.sqrt(varX1 + 4.0*avgEnStageSimpleInf*avgEnStageSimpleInf*varX0 - 4*avgEnStageSimpleInf*dataErrStageSimpleInf.getValue(0)*dataErrStageSimpleInf.getValue(1)*corX0X1);

        //8 HMA EC stage
        DataGroup dataStageECInf = (DataGroup) accumulatorStageECInf.getData();
        IData dataAvgStageECInf = dataStageECInf.getData(accumulatorStageECInf.AVERAGE.index);
        IData dataErrStageECInf = dataStageECInf.getData(accumulatorStageECInf.ERROR.index);
        IData dataCorStageECInf = dataStageECInf.getData(accumulatorStageECInf.BLOCK_CORRELATION.index);
        IData dataCovStageECInf = dataStageECInf.getData(accumulatorStageECInf.COVARIANCE.index);
        double avgEnStageECInf = dataAvgStageECInf.getValue(0);
        double errEnStageECInf = dataErrStageECInf.getValue(0);
        double corEnStageECInf = dataCorStageECInf.getValue(0);
        System.out.println(" En_stage_ec:      " + avgEnStageECInf + "   err: " + errEnStageECInf + " cor: " + corEnStageECInf);
        double CvnStageEC  = kB_beta2*(dataAvgStageECInf.getValue(1) - avgEnStageECInf*avgEnStageECInf);
        varX0 = errEnStageECInf*errEnStageECInf;
        varX1 = dataErrStageECInf.getValue(1)*dataErrStageECInf.getValue(1);
        corX0X1 = dataCovStageECInf.getValue(1)/Math.sqrt(dataCovStageECInf.getValue(0))/Math.sqrt(dataCovStageECInf.getValue(3));
        double errCvnStageEC = kB_beta2*Math.sqrt(varX1 + 4.0*avgEnStageECInf*avgEnStageECInf*varX0 - 4*avgEnStageECInf*errEnStageECInf*dataErrStageECInf.getValue(1)*corX0X1);
        if (errEnStageECInf < 1e-10){
            errCvnStageEC = kB_beta2*Math.sqrt(varX1 + 4.0*avgEnStageECInf*avgEnStageECInf*varX0);
        }

        System.out.println();
        System.out.println(" Cvn_prim:         " + CvnPrimInf +          "   err: " + errCvnPrimInf);
        System.out.println(" Cvn_vir:          " + CvnVirInf +           "   err: " + errCvnVirInf);
        System.out.println(" Cvn_cvir:         " + CvnCentVirInf +       "   err: " + errCvnCentVirInf);
        System.out.println(" Cvn_hmac:         " + CvnHMAcInf +          "   err: " + errCvnHMAcInf);
        System.out.println(" Cvn_nm_simple:    " + Cvn_nm_simpleInf +    "   err: " + errCvnNMsimpleInf);
        System.out.println(" Cvn_nm_ec:        " + CvnNMECInf +          "   err: " + errCvnNMECInf);
        System.out.println(" Cvn_stage_simple: " + Cvn_stage_simple + "   err: " + errCvnStageSimple);
        System.out.println(" Cvn_stage_ec:     " + CvnStageEC +       "   err: " + errCvnStageEC);


        long endTime = System.currentTimeMillis();
        System.out.println("\n time: (min) " + (endTime - startTime)/60.0/1000.0);
    }

    public enum MoveChoice {Real, NM, NMEC, Stage, StageEC};

    public static class OctaneParams extends ParameterBase {
        public double temperature = 1;
        public double hbar = 1;
        public double omega = 1;
        public double k4 = 0;
        public long steps = 10_000_000;
        public boolean isGraphic = false;
        public boolean isTIA = false;
        public boolean zerok0 = false;
        public boolean onlyCentroid = false;
        public double mass = 1.0;
        public MoveChoice coordType = MoveChoice.Real;
        public double timeStep = -1;
        public int nBeads = -1;
        public int nShifts = 0;

    }
}