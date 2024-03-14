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
        moveStageEC = new MCMoveHOReal2(space, pmAgg, random, temperature, omega*omega, box, hbar);
        moveStageSimple = new MCMoveHOReal2(space, pmAgg, random, temperature, 0, box, hbar);

        if (coordType == MoveChoice.Real) {
            integrator = new IntegratorLangevin(pmAgg, random, timeStep, temperature, box, gammaLangevin);
        } else if (coordType == MoveChoice.NM) {
            integrator = new IntegratorLangevinPINMInf(pmAgg, random, timeStep, temperature, box, gammaLangevin, hbar, omega, 0);
        } else if (coordType == MoveChoice.NMEC) {
            integrator = new IntegratorLangevinPINMInf(pmAgg, random, timeStep, temperature, box, gammaLangevin, hbar, omega, omegaSample);
        } else if (coordType == MoveChoice.Stage) {
            integrator = new IntegratorLangevinPIInf(pmAgg, random, timeStep, temperature, box, gammaLangevin, hbar, omega, 0, moveStageEC);
        } else { //StageEC -- default
            integrator = new IntegratorLangevinPIInf(pmAgg, random, timeStep, temperature, box, gammaLangevin,  hbar, omega, omegaSample, moveStageEC);
        }

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
            params.steps = 1000000;
            params.temperature = 0.1;
            params.omega = 1;
            params.k4 = 0.1;
            params.coordType = MoveChoice.Real;
//            params.coordType = MoveChoice.NMEC;
//            params.coordType = MoveChoice.StageEC;
            params.onlyCentroid = true;
            params.timeStep = 0.01;
        }

        int nShifts = params.nShifts;
        double mass = params.mass;
        int nBeads = params.nBeads;
        double temperature = params.temperature;
        double beta = 1 / temperature;
        double hbar = params.hbar;
        double omega = params.omega;
        double k4 = params.k4;
        boolean isGraphic = params.isGraphic;
        long steps = params.steps;
        long stepsEq = steps / 10;
        boolean isTIA = params.isTIA;
        boolean zerok0 = params.zerok0;
        boolean onlyCentroid = params.onlyCentroid;
        MoveChoice coordType = params.coordType;
        double omega2 = omega * omega;

        double x = hbar * omega / temperature;
        if (nBeads == -1) {
            nBeads = (int) (20 * x); //20*x and 30*x are good for HO and AO, resp.
        }
        double alpha = beta * hbar * omega / nBeads;
        double mOmegaF2 = nBeads == 1 ? 0 : mass * omega2 / nBeads / alpha / Math.sinh(alpha);
        double mOmegaH2 = 2 * mass * omega2 * Math.tanh(alpha / 2) / nBeads / alpha;
        double omegaRing = Math.sqrt(mOmegaH2 * nBeads / mass);
        double omegaBead = Math.sqrt(mOmegaF2 * nBeads / mass);

        double omegaN = Math.sqrt(nBeads) * temperature / hbar;
        double timeStep = params.timeStep;
        if (timeStep == -1) {
            double c = 0.1;
            if (params.coordType == MoveChoice.Real) {
                timeStep = c / omegaBead;// mi=m/n in real space, so m wn^2 = (m/n)*(n wn^2)==> dt~1/[wn sqrt(n)]
            } else if (params.coordType == MoveChoice.NM) {
                double s = omega2 / nBeads;
                timeStep = c * Math.sqrt(s) / omega; // which is 1/sqrt(n)
            } else if (params.coordType == MoveChoice.Stage) {
                double s = omega2 / omegaN / omegaN * (1 + 1.0 / 12.0 * (nBeads * nBeads - 1.0) / nBeads);
                timeStep = c * Math.sqrt(s) / omega; // for large n, timeStep ~ hbar/T
            } else {
                timeStep = c / omegaRing;
            }
        }

//        if (isTIA){
//            omega2 = omega2*(1.0 + omega2/12.0/(nBeads*omegaN2));
//        }
//        if (zerok0) omega2 = 0;

        double gammaLangevin = 2 * omegaRing;
        final SimQuantumAOPIMDInf sim = new SimQuantumAOPIMDInf(Space1D.getInstance(), coordType, mass, timeStep, gammaLangevin, nBeads, temperature, k4, omega, mOmegaF2, mOmegaH2, isTIA, hbar);
        sim.integrator.reset();
        System.out.println(" PIMDinf-" + coordType);
        System.out.println(" mass: " + mass);
        System.out.println(" T: " + temperature);
        System.out.println(" hbar: " + hbar);
        System.out.println(" w: " + omega);
        System.out.println(" mOmegaF2: " + mOmegaF2 + " , mOmegaH2: " + mOmegaH2);
        System.out.println(" x = beta*hbar*w = " + hbar * omega / temperature);
        System.out.println(" nBeads: " + nBeads);
        System.out.println(" nShifts: " + nShifts);
        System.out.println(" steps: " + steps + " stepsEq: " + stepsEq);
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
        if (isTIA) {
//            meterPrim = new MeterPIPrim(sim.pmBonding, sim.pcP1EnTIA, nBeads, sim.betaN);
//            meterVir = new MeterPIVirTIA(sim.pcP1EnTIA, sim.pcP1, sim.betaN, nBeads, sim.box);
//            meterCentVir = new MeterPICentVirTIA(sim.pcP1EnTIA, sim.pcP1, sim.betaN, nBeads, sim.box);
            meterHMAcInf = null;
//            meterHMA = new MeterPIHMATIA(sim.pmBonding, sim.pcP1EnTIA, sim.pcP1, sim.betaN, nBeads, omega2, sim.box);
//            meterHMAcent = null;
        } else {
            meterCentVirInf = new MeterPICentVirInf(sim.pcP1harm, sim.pcP1ah, temperature, sim.box, nBeads, omega, hbar);
            if (!onlyCentroid) {
                meterPrimInf = new MeterPIPrimInf(sim.pmBonding, sim.pcP1harm, sim.pcP1ah, nBeads, temperature, sim.box, omega, hbar);
                meterVirInf = new MeterPIVirInf(sim.pcP1harm, sim.pcP1ah, temperature, sim.box, nBeads, omega, hbar);
                meterVirHarmStagingInf = new MeterPIVirHarmStagingInf(sim.pcP1ah, temperature, sim.box, nBeads, omega, hbar);
                meterHMAcInf = new MeterPIHMAcInf(sim.pcP1harm, sim.pcP1ah, temperature, sim.box, nBeads, omega, hbar);
                meterNMSimpleInf = new MeterPIHMAInf(sim.pmBonding, sim.pcP1harm, sim.pcP1ah, temperature, nBeads, omega, 0, sim.box, hbar);
                meterNMECInf = new MeterPIHMAInf(sim.pmBonding, sim.pcP1harm, sim.pcP1ah, temperature, nBeads, omega, omega, sim.box, hbar);
                meterStageSimpleInf = new MeterPIHMAReal2Inf(sim.pmBonding, sim.pcP1harm, sim.pcP1ah, temperature, nBeads, omega, 0, sim.box, hbar);
                meterStageSimpleInf.setNumShifts(nShifts);
                meterStageECInf = new MeterPIHMAReal2Inf(sim.pmBonding, sim.pcP1harm, sim.pcP1ah, temperature, nBeads, omega, omega, sim.box, hbar);
                meterStageECInf.setNumShifts(nShifts);
            }
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
            simGraphic.makeAndDisplayFrame("PIMD - " + coordType);
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
        if (meterPrimInf != null) {
            DataPumpListener accumulatorPumpPrimInf = new DataPumpListener(meterPrimInf, accumulatorPrimInf, interval);
            sim.integrator.getEventManager().addListener(accumulatorPumpPrimInf);
        }
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
        System.out.println("\n T_measured: " + avgKE / (0.5 * sim.dim * nBeads) + " +/- " + errKE / (0.5 * sim.dim * nBeads) + " cor: " + corKE);

        double kB_beta2 = beta * beta;
        double varX0, varX1, corX0X1;


        // Cent-Vir
        DataGroup dataCentVir = (DataGroup) accumulatorCentVirInf.getData();
        IData dataAvgCentVir = dataCentVir.getData(accumulatorCentVirInf.AVERAGE.index);
        IData dataErrCentVir = dataCentVir.getData(accumulatorCentVirInf.ERROR.index);
        IData dataCorCentVir = dataCentVir.getData(accumulatorCentVirInf.BLOCK_CORRELATION.index);
        IData dataCovCentVir = dataCentVir.getData(accumulatorCentVirInf.COVARIANCE.index);
        double avgEnCentVir = dataAvgCentVir.getValue(0);
        double errEnCentVir = dataErrCentVir.getValue(0);
        double corEnCentVir = dataCorCentVir.getValue(0);
//        double CvnCentVir = kB_beta2*(dataAvgCentVir.getValue(1) - avgEnCentVir*avgEnCentVir);
//        varX0 = errEnCentVir*errEnCentVir;
//        varX1 = dataErrCentVir.getValue(1)*dataErrCentVir.getValue(1);
//        corX0X1 = dataCovCentVir.getValue(1)/Math.sqrt(dataCovCentVir.getValue(0))/Math.sqrt(dataCovCentVir.getValue(3));
//        double errCvnCentVir = kB_beta2*Math.sqrt(varX1 + 4.0*avgEnCentVir*avgEnCentVir*varX0 - 4*avgEnCentVir*dataErrCentVir.getValue(0)*dataErrCentVir.getValue(1)*corX0X1);


        double avgEnPrim = 0, avgEnVir = 0, avgEnHMAc = 0, avgEnNMSimple = 0, avgEnNMEC = 0, avgEnStageSimple = 0, avgEnStageEC = 0;
        double errEnPrim = 0, errEnVir = 0, errEnHMAc = 0, errEnNMSimple = 0, errEnNMEC = 0, errEnStageSimple = 0, errEnStageEC = 0;
        double corEnPrim = 0, corEnVir = 0, corEnHMAc = 0, corEnNMSimple = 0, corEnNMEC = 0, corEnStageSimple = 0, corEnStageEC = 0;
        double CvnPrim = 0, CvnVir = 0, CvnHMAc = 0, CvnNMSimple = 0, CvnNMEC = 0, CvnStageSimple = 0, CvnStageEC = 0;
        double errCvnPrim = 0, errCvnVir = 0, errCvnHMAc = 0, errCvnNMsimple = 0, errCvnNMEC = 0, errCvnStageSimple = 0, errCvnStageEC = 0;

        if (!onlyCentroid) {
            // Prim
            DataGroup dataPrim = (DataGroup) accumulatorPrimInf.getData();
            IData dataAvgPrim = dataPrim.getData(accumulatorPrimInf.AVERAGE.index);
            IData dataErrPrim = dataPrim.getData(accumulatorPrimInf.ERROR.index);
            IData dataCorPrim = dataPrim.getData(accumulatorPrimInf.BLOCK_CORRELATION.index);
            IData dataCovPrim = dataPrim.getData(accumulatorPrimInf.COVARIANCE.index);
            avgEnPrim = dataAvgPrim.getValue(0);
            errEnPrim = dataErrPrim.getValue(0);
            corEnPrim = dataCorPrim.getValue(0);
            CvnPrim = kB_beta2 * (dataAvgPrim.getValue(1) - avgEnPrim * avgEnPrim);
            varX0 = errEnPrim * errEnPrim;
            varX1 = dataErrPrim.getValue(1) * dataErrPrim.getValue(1);
            corX0X1 = dataCovPrim.getValue(1) / Math.sqrt(dataCovPrim.getValue(0)) / Math.sqrt(dataCovPrim.getValue(3));
            errCvnPrim = kB_beta2 * Math.sqrt(varX1 + 4.0 * avgEnPrim * avgEnPrim * varX0 - 4 * avgEnPrim * dataErrPrim.getValue(0) * dataErrPrim.getValue(1) * corX0X1);

            // Vir
            DataGroup dataVir = (DataGroup) accumulatorVirInf.getData();
            IData dataAvgVir = dataVir.getData(accumulatorVirInf.AVERAGE.index);
            IData dataErrVir = dataVir.getData(accumulatorVirInf.ERROR.index);
            IData dataCorVir = dataVir.getData(accumulatorVirInf.BLOCK_CORRELATION.index);
            IData dataCovVir = dataVir.getData(accumulatorVirInf.COVARIANCE.index);
            avgEnVir = dataAvgVir.getValue(0);
            errEnVir = dataErrVir.getValue(0);
            corEnVir = dataCorVir.getValue(0);
            CvnVir = kB_beta2 * (dataAvgVir.getValue(1) - avgEnVir * avgEnVir);
            varX0 = errEnVir * errEnVir;
            varX1 = dataErrVir.getValue(1) * dataErrVir.getValue(1);
            corX0X1 = dataCovVir.getValue(1) / Math.sqrt(dataCovVir.getValue(0)) / Math.sqrt(dataCovVir.getValue(3));
            errCvnVir = kB_beta2 * Math.sqrt(varX1 + 4.0 * avgEnVir * avgEnVir * varX0 - 4 * avgEnVir * dataErrVir.getValue(0) * dataErrVir.getValue(1) * corX0X1);

            // HMAc
            DataGroup dataHMAc = (DataGroup) accumulatorHMAcInf.getData();
            IData dataAvgHMAc = dataHMAc.getData(accumulatorHMAcInf.AVERAGE.index);
            IData dataErrHMAc = dataHMAc.getData(accumulatorHMAcInf.ERROR.index);
            IData dataCorHMAc = dataHMAc.getData(accumulatorHMAcInf.BLOCK_CORRELATION.index);
            IData dataCovHMAc = dataHMAc.getData(accumulatorHMAcInf.COVARIANCE.index);
            avgEnHMAc = dataAvgHMAc.getValue(0);
            errEnHMAc = dataErrHMAc.getValue(0);
            corEnHMAc = dataCorHMAc.getValue(0);
            CvnHMAc = kB_beta2 * (dataAvgHMAc.getValue(1) - avgEnHMAc * avgEnHMAc);
            varX0 = errEnHMAc * errEnHMAc;
            varX1 = dataErrHMAc.getValue(1) * dataErrHMAc.getValue(1);
            corX0X1 = dataCovHMAc.getValue(1) / Math.sqrt(dataCovHMAc.getValue(0)) / Math.sqrt(dataCovHMAc.getValue(3));
            errCvnHMAc = kB_beta2 * Math.sqrt(varX1 + 4.0 * avgEnHMAc * avgEnHMAc * varX0 - 4 * avgEnHMAc * dataErrHMAc.getValue(0) * dataErrHMAc.getValue(1) * corX0X1);
            if (errEnHMAc < 1e-10) {
                errCvnHMAc = kB_beta2 * Math.sqrt(varX1 + 4.0 * avgEnHMAc * avgEnHMAc * varX0);
            }


            // HMA simple NM
            DataGroup dataNMsimple = (DataGroup) accumulatorNMSimpleInf.getData();
            IData dataAvgNMsimple = dataNMsimple.getData(accumulatorNMSimpleInf.AVERAGE.index);
            IData dataErrNMsimple = dataNMsimple.getData(accumulatorNMSimpleInf.ERROR.index);
            IData dataCorNMsimple = dataNMsimple.getData(accumulatorNMSimpleInf.BLOCK_CORRELATION.index);
            IData dataCovNMsimple = dataNMsimple.getData(accumulatorNMSimpleInf.COVARIANCE.index);
            avgEnNMSimple = dataAvgNMsimple.getValue(0);
            errEnNMSimple = dataErrNMsimple.getValue(0);
            corEnNMSimple = dataCorNMsimple.getValue(0);
            CvnNMSimple = kB_beta2 * (dataAvgNMsimple.getValue(1) - avgEnNMSimple * avgEnNMSimple);
            varX0 = errEnNMSimple * errEnNMSimple;
            varX1 = dataErrNMsimple.getValue(1) * dataErrNMsimple.getValue(1);
            corX0X1 = dataCovNMsimple.getValue(1) / Math.sqrt(dataCovNMsimple.getValue(0)) / Math.sqrt(dataCovNMsimple.getValue(3));
            errCvnNMsimple = kB_beta2 * Math.sqrt(varX1 + 4.0 * avgEnNMSimple * avgEnNMSimple * varX0 - 4 * avgEnNMSimple * dataErrNMsimple.getValue(0) * dataErrNMsimple.getValue(1) * corX0X1);

            // HMA EC NM
            DataGroup dataNMEC = (DataGroup) accumulatorNMECInf.getData();
            IData dataAvgNMEC = dataNMEC.getData(accumulatorNMECInf.AVERAGE.index);
            IData dataErrNMEC = dataNMEC.getData(accumulatorNMECInf.ERROR.index);
            IData dataCorNMEC = dataNMEC.getData(accumulatorNMECInf.BLOCK_CORRELATION.index);
            IData dataCovNMEC = dataNMEC.getData(accumulatorNMECInf.COVARIANCE.index);
            avgEnNMEC = dataAvgNMEC.getValue(0);
            errEnNMEC = dataErrNMEC.getValue(0);
            corEnNMEC = dataCorNMEC.getValue(0);
            CvnNMEC = kB_beta2 * (dataAvgNMEC.getValue(1) - avgEnNMEC * avgEnNMEC);
            varX0 = errEnNMEC * errEnNMEC;
            varX1 = dataErrNMEC.getValue(1) * dataErrNMEC.getValue(1);
            corX0X1 = dataCovNMEC.getValue(1) / Math.sqrt(dataCovNMEC.getValue(0)) / Math.sqrt(dataCovNMEC.getValue(3));
            errCvnNMEC = kB_beta2 * Math.sqrt(varX1 + 4.0 * avgEnNMEC * avgEnNMEC * varX0 - 4 * avgEnNMEC * errEnNMEC * dataErrNMEC.getValue(1) * corX0X1);
            if (errEnNMEC < 1e-10) {
                errCvnNMEC = kB_beta2 * Math.sqrt(varX1 + 4.0 * avgEnNMEC * avgEnNMEC * varX0);
            }

            // HMA simple stage
            DataGroup dataStageSimple = (DataGroup) accumulatorStageSimpleInf.getData();
            IData dataAvgStageSimple = dataStageSimple.getData(accumulatorStageSimpleInf.AVERAGE.index);
            IData dataErrStageSimple = dataStageSimple.getData(accumulatorStageSimpleInf.ERROR.index);
            IData dataCorStageSimple = dataStageSimple.getData(accumulatorStageSimpleInf.BLOCK_CORRELATION.index);
            IData dataCovStageSimple = dataStageSimple.getData(accumulatorStageSimpleInf.COVARIANCE.index);
            avgEnStageSimple = dataAvgStageSimple.getValue(0);
            errEnStageSimple = dataErrStageSimple.getValue(0);
            corEnStageSimple = dataCorStageSimple.getValue(0);
            CvnStageSimple = kB_beta2 * (dataAvgStageSimple.getValue(1) - avgEnStageSimple * avgEnStageSimple);
            varX0 = errEnStageSimple * errEnStageSimple;
            varX1 = dataErrStageSimple.getValue(1) * dataErrStageSimple.getValue(1);
            corX0X1 = dataCovStageSimple.getValue(1) / Math.sqrt(dataCovStageSimple.getValue(0)) / Math.sqrt(dataCovStageSimple.getValue(3));
            errCvnStageSimple = kB_beta2 * Math.sqrt(varX1 + 4.0 * avgEnStageSimple * avgEnStageSimple * varX0 - 4 * avgEnStageSimple * dataErrStageSimple.getValue(0) * dataErrStageSimple.getValue(1) * corX0X1);

            // HMA EC stage
            DataGroup dataStageEC = (DataGroup) accumulatorStageECInf.getData();
            IData dataAvgStageEC = dataStageEC.getData(accumulatorStageECInf.AVERAGE.index);
            IData dataErrStageEC = dataStageEC.getData(accumulatorStageECInf.ERROR.index);
            IData dataCorStageEC = dataStageEC.getData(accumulatorStageECInf.BLOCK_CORRELATION.index);
            IData dataCovStageEC = dataStageEC.getData(accumulatorStageECInf.COVARIANCE.index);
            avgEnStageEC = dataAvgStageEC.getValue(0);
            errEnStageEC = dataErrStageEC.getValue(0);
            corEnStageEC = dataCorStageEC.getValue(0);
            CvnStageEC = kB_beta2 * (dataAvgStageEC.getValue(1) - avgEnStageEC * avgEnStageEC);
            varX0 = errEnStageEC * errEnStageEC;
            varX1 = dataErrStageEC.getValue(1) * dataErrStageEC.getValue(1);
            corX0X1 = dataCovStageEC.getValue(1) / Math.sqrt(dataCovStageEC.getValue(0)) / Math.sqrt(dataCovStageEC.getValue(3));
            errCvnStageEC = kB_beta2 * Math.sqrt(varX1 + 4.0 * avgEnStageEC * avgEnStageEC * varX0 - 4 * avgEnStageEC * errEnStageEC * dataErrStageEC.getValue(1) * corX0X1);
            if (errEnStageEC < 1e-10) {
                errCvnStageEC = kB_beta2 * Math.sqrt(varX1 + 4.0 * avgEnStageEC * avgEnStageEC * varX0);
            }
        }

        if (onlyCentroid) {
            System.out.println(" En_cvir:          " + avgEnCentVir + "   err: " + errEnCentVir + " cor: " + corEnCentVir);
            System.out.println();
//            System.out.println(" Cvn_prim:         " + CvnPrim +          "   err: " + errCvnPrim);
//            System.out.println(" Cvn_cvir:         " + CvnCentVir +       "   err: " + errCvnCentVir);
        } else {
            System.out.println("\n En_prim:          " + avgEnPrim + "   err: " + errEnPrim + " cor: " + corEnPrim);
            System.out.println(" En_vir:           " + avgEnVir + "   err: " + errEnVir + " cor: " + corEnVir);
            System.out.println(" En_cvir:          " + avgEnCentVir + "   err: " + errEnCentVir + " cor: " + corEnCentVir);
            System.out.println(" En_hmac:          " + avgEnHMAc + "   err: " + errEnHMAc + " cor: " + corEnHMAc);
            System.out.println(" En_nm_simple:     " + avgEnNMSimple + "   err: " + errEnNMSimple + " cor: " + corEnNMSimple);
            System.out.println(" En_nm_ec:         " + avgEnNMEC + "   err: " + errEnNMEC + " cor: " + corEnNMEC);
            System.out.println(" En_stage_simple:  " + avgEnStageSimple + "   err: " + errEnStageSimple + " cor: " + corEnStageSimple);
            System.out.println(" En_stage_ec:      " + avgEnStageEC + "   err: " + errEnStageEC + " cor: " + corEnStageEC);
            System.out.println();
            System.out.println(" Cvn_prim:         " + CvnPrim + "   err: " + errCvnPrim);
            System.out.println(" Cvn_vir:          " + CvnVir + "   err: " + errCvnVir);
//            System.out.println(" Cvn_cvir:         " + CvnCentVir +       "   err: " + errCvnCentVir);
            System.out.println(" Cvn_hmac:         " + CvnHMAc + "   err: " + errCvnHMAc);
            System.out.println(" Cvn_nm_simple:    " + CvnNMSimple + "   err: " + errCvnNMsimple);
            System.out.println(" Cvn_nm_ec:        " + CvnNMEC + "   err: " + errCvnNMEC);
            System.out.println(" Cvn_stage_simple: " + CvnStageSimple + "   err: " + errCvnStageSimple);
            System.out.println(" Cvn_stage_ec:     " + CvnStageEC + "   err: " + errCvnStageEC);
        }



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
        public boolean onlyCentroid = true;
        public double mass = 1.0;
        public MoveChoice coordType = MoveChoice.Real;
        public double timeStep = -1;
        public int nBeads = -1;
        public int nShifts = 0;

    }
}