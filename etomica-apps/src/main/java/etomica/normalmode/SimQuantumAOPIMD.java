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

public class SimQuantumAOPIMD extends Simulation {

    public PotentialComputeField pcP1, pcP1EnTIA;
    public final IntegratorMD integrator;
    public final PotentialMasterBonding pmBonding;
    public final Box box;
    public final PotentialComputeAggregate pmAgg;
    public IPotential1 p1ah;
    public P1AnharmonicTIA p1ahUeff;
    public double betaN;
    public double k2_kin;
    public int nBeads;
    public MCMoveHOReal2 moveStageSimple, moveStageEC;
    public int dim;

    public SimQuantumAOPIMD(Space space, MoveChoice coordType, double mass, double timeStep, double gammaLangevin, int nBeads, double temperature, double k2, double k4, double omega2, boolean isTIA, double hbar) {
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
        double beta = 1.0/temperature;
        betaN = beta/nBeads;
        double omegaN = Math.sqrt(nBeads)/(hbar*beta);
        k2_kin = nBeads == 1 ? 0 : mass*omegaN*omegaN;

        // integrate analytically to get infinite-n
//        if (false){
//            double w = Math.sqrt(omega2);
//            double R = betaN*hbar*w;
//            k2_kin = nBeads == 1 ? 0 : mass*omegaN*omegaN*R/Math.tanh(R);
//            k2 = 2*nBeads*mass*omegaN*omegaN*R*Math.tanh(R/2);
//        }

        P2Harmonic p2Bond = new P2Harmonic(k2_kin, 0);
        List<int[]> pairs = new ArrayList<>();
        for (int i=0; i<nBeads; i++) {
            int[] p = new int[]{i,(i+1)%nBeads};
            pairs.add(p);
        }
        pmBonding.setBondingPotentialPair(species, p2Bond, pairs);

        pcP1 = new PotentialComputeField(getSpeciesManager(), box);

        if (isTIA){
            double facUeff = 1.0;
            p1ahUeff = new P1AnharmonicTIA(space, k2, k4, nBeads, mass*omegaN*omegaN, facUeff);
            pcP1.setFieldPotential(species.getLeafType(), p1ahUeff);
        } else {
            p1ah = new P1Anharmonic(space, k2/nBeads, k4/nBeads);
            pcP1.setFieldPotential(species.getLeafType(), p1ah);
        }

        PotentialComputeAggregate.localStorageDefault = true;
        pmAgg = new PotentialComputeAggregate(pmBonding, pcP1);

        double facEn = 3.0;
        P1AnharmonicTIA p1ahEn = new P1AnharmonicTIA(space, k2, k4, nBeads, mass*omegaN*omegaN, facEn);
        pcP1EnTIA = new PotentialComputeField(getSpeciesManager(), box);
        pcP1EnTIA.setFieldPotential(species.getLeafType(), p1ahEn);

        if (coordType == MoveChoice.Real) {
            integrator = new IntegratorLangevin(pmAgg, random, timeStep, temperature, box, gammaLangevin);
        } else if (coordType == MoveChoice.NM) {
            MCMoveHO move = new MCMoveHO(space, pmAgg, random, temperature, 0, box, hbar);
            integrator = new IntegratorLangevinPINM(pmAgg, random, timeStep, temperature, box, gammaLangevin, move, hbar, omega2);
        } else if (coordType == MoveChoice.NMEC) {
            MCMoveHO move = new MCMoveHO(space, pmAgg, random, temperature, omega2, box, hbar);
            integrator = new IntegratorLangevinPINM(pmAgg, random, timeStep, temperature, box, gammaLangevin, move, hbar, omega2);
        } else if (coordType == MoveChoice.Stage) {
            MCMoveHOReal2 move = new MCMoveHOReal2(space, pmAgg, random, temperature, 0, box, hbar);
            integrator = new IntegratorLangevinPI(pmAgg, random, timeStep, temperature, box, gammaLangevin, move, hbar, omega2);
        } else { //StageEC -- default
            MCMoveHOReal2 move = new MCMoveHOReal2(space, pmAgg, random, temperature, omega2, box, hbar);
            integrator = new IntegratorLangevinPI(pmAgg, random, timeStep, temperature, box, gammaLangevin, move, hbar, omega2);
        }

        moveStageEC = new MCMoveHOReal2(space, pmAgg, random, temperature, omega2, box, hbar);
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
            params.steps = 1000000;
            params.temperature = 1;
            params.k2 = 1;
            params.k4 = 24;
//            params.coordType = MoveChoice.Real;
//            params.coordType = MoveChoice.NM;
            params.coordType = MoveChoice.Stage;
            params.nBeads = 10;
            params.nShifts = 0;


            params.timeStep = 0.001;



        }

        int nShifts = params.nShifts;
        double mass = params.mass;
        double temperature = params.temperature;
        double hbar = params.hbar;
        double k2 = params.k2;
        double k4 = params.k4;
        double gammaLangevin = 2*Math.sqrt(k2);
        boolean isGraphic = params.isGraphic;
        long steps = params.steps;
        long stepsEq = steps/10;
        boolean isTIA = params.isTIA;
        boolean zerok0 = params.zerok0;
        boolean onlyCentroid = params.onlyCentroid;
        MoveChoice coordType = params.coordType;
        double omega = Math.sqrt(k2/mass);
        double omega2 = k2/mass;

        double x = 1/temperature*hbar*omega;
        int nBeads = params.nBeads;
        if (nBeads == -1){
            nBeads = (int) (30*x); //20*x and 30*x are good for HO and AO, resp.
        }

        double omegaN = Math.sqrt(nBeads)*temperature/hbar;
        double timeStep = params.timeStep;
        if (timeStep == -1) {
            double c = 0.1;
            if (params.coordType == MoveChoice.Real) {
                timeStep = c/omegaN/Math.sqrt(nBeads);//c/(20*w)
            } else if (params.coordType == MoveChoice.Stage) {
                double s = 1.0 + 1.0/12.0*omega2/omegaN/omegaN/nBeads*(nBeads*nBeads-1);
                timeStep = c*Math.sqrt(s)/omega;
            } else {
                timeStep = c/omega;
            }
        }

        if (isTIA){
            omega2 = omega2*(1.0 + omega2/12.0/(nBeads*omegaN*omegaN));
        }
//        if (zerok0) omega2 = 0;

        final SimQuantumAOPIMD sim = new SimQuantumAOPIMD(Space1D.getInstance(), coordType, mass, timeStep, gammaLangevin, nBeads, temperature, k2, k4, omega2, isTIA, hbar);
        sim.integrator.reset();

        System.out.println(" PIMD-" + coordType);
        System.out.println(" mass: " + mass);
        System.out.println(" T: " + temperature);
        System.out.println(" hbar: " + hbar);
        System.out.println(" w: " + Math.sqrt(omega2));
        System.out.println(" wn: " + omegaN  + " , w/sqrt(n): " + Math.sqrt(omega2)/Math.sqrt(nBeads));
        System.out.println(" x = beta*hbar*w = " + hbar*omega/temperature);
        System.out.println(" nBeads: " + nBeads);
        System.out.println(" nShifts: "+ nShifts);
        System.out.println(" timestep: " + timeStep);

        System.out.println(" steps: " +  steps + " stepsEq: " + stepsEq);
        System.out.println(" k2: " + k2);
        System.out.println(" k4: " + k4);
        System.out.println(" isTIA: " + isTIA);
        System.out.println(" gammaLangevin: " + gammaLangevin);

        System.out.println("\n Quantum Harmonic Oscillator Theory");
        System.out.println(" ====================================");
        double alpha = 1 + 0.5*Math.pow(hbar*sim.betaN*omega,2)+0.5*hbar* sim.betaN*omega*Math.sqrt(4+Math.pow(hbar* sim.betaN*omega,2));
        double alpha2 = alpha*alpha;
        double hbar2 = hbar*hbar;
        double dAlphaDBeta = 2.0/temperature*hbar2*omega2/nBeads/nBeads*alpha2/(alpha2-1);
        double dAlphadT = -1.0/temperature/temperature*dAlphaDBeta;
        double EnQ = sim.space.D()*(hbar2*omega2)*sim.betaN*alpha/(alpha*alpha-1)*(Math.pow(alpha,nBeads)+1)/(Math.pow(alpha,nBeads)-1);
        double numerator = 1 + alpha2 - Math.pow(alpha,2*nBeads)*(alpha2+1)-2*nBeads*(alpha2-1)*Math.pow(alpha,nBeads);
        double denominator = (alpha2-1)*(alpha2-1)*(Math.pow(alpha,nBeads)-1)*(Math.pow(alpha,nBeads)-1);
        double CvnQ = sim.space.D()*hbar2*omega2*sim.betaN*dAlphadT*numerator/denominator-1/temperature/temperature*EnQ*temperature;
        double EnQinf = sim.space.D()*hbar*omega*(0.5 + 1/(Math.exp(nBeads*sim.betaN*hbar*omega)-1.0));
        double CvnQinf = sim.space.D()*Math.pow(1.0/temperature*hbar*omega/2/Math.sinh(1.0/temperature*hbar*omega/2), 2);
        double EnC = sim.space.D()*temperature;
        double CvnC = sim.space.D();
        System.out.println(" En_ho_c: " + EnC);
        System.out.println(" En_ho_q: " + EnQ);
        System.out.println(" E_ho_q: " + EnQinf);
        System.out.println(" Cvn_ho_c: " + CvnC);
        System.out.println(" Cvn_ho_q: " + CvnQ);
        System.out.println(" Cv_ho_c: " + CvnQinf + "\n");

        DataSourceScalar meterKE = sim.integrator.getMeterKineticEnergy();

//        MeterMSDHO meterMSDHO = new MeterMSDHO(sim.box);
        MeterPIPrim meterPrim = null;
        MeterPIVir meterVir = null;
        MeterPICentVir meterCentVir = null;
        MeterPIHMAc meterHMAc = null;
        MeterPIHMA meterNMEC = null;
        MeterPIHMA meterNMSimple = null;
        MeterPIHMAReal2 meterStageEC = null;
        MeterPIHMAReal2 meterStageSimple = null;
        if (isTIA){
//            meterPrim = new MeterPIPrim(sim.pmBonding, sim.pcP1EnTIA, nBeads, sim.betaN);
//            meterVir = new MeterPIVirTIA(sim.pcP1EnTIA, sim.pcP1, sim.betaN, nBeads, sim.box);
//            meterCentVir = new MeterPICentVirTIA(sim.pcP1EnTIA, sim.pcP1, sim.betaN, nBeads, sim.box);
            meterHMAc = null;
//            meterHMA = new MeterPIHMATIA(sim.pmBonding, sim.pcP1EnTIA, sim.pcP1, sim.betaN, nBeads, omega2, sim.box);
//            meterHMAcent = null;
        } else {
            meterPrim = new MeterPIPrim(sim.pmBonding, sim.pcP1, nBeads, temperature, sim.box);
            meterVir = new MeterPIVir(sim.pcP1, temperature, sim.box);
            meterCentVir = new MeterPICentVir(sim.pcP1, temperature, nBeads, sim.box);
            meterHMAc = new MeterPIHMAc(sim.pcP1, temperature, nBeads, sim.box);
            meterNMSimple = new MeterPIHMA(sim.pmBonding, sim.pcP1, sim.betaN, nBeads, 0, sim.box, hbar);
            meterNMEC = new MeterPIHMA(sim.pmBonding, sim.pcP1, sim.betaN, nBeads, omega2, sim.box, hbar);
            meterStageSimple = new MeterPIHMAReal2(sim.pmBonding, sim.pcP1, nBeads, temperature, sim.moveStageSimple);
            meterStageSimple.setNumShifts(nShifts);
            meterStageEC = new MeterPIHMAReal2(sim.pmBonding, sim.pcP1, nBeads, temperature, sim.moveStageEC);
            meterStageEC.setNumShifts(nShifts);
        }

        MeterPIHMAvir meterHMAvir = new MeterPIHMAvir(sim.pmBonding, sim.pcP1, sim.betaN, nBeads, omega2, sim.box, hbar);//Bad!!

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
        AccumulatorAverageCovariance accumulatorPrim = new AccumulatorAverageCovariance(blockSize);
        if (meterPrim != null) {
            DataPumpListener accumulatorPumpPrim = new DataPumpListener(meterPrim, accumulatorPrim, interval);
            sim.integrator.getEventManager().addListener(accumulatorPumpPrim);
        }

        //2 Virial
        AccumulatorAverageCovariance accumulatorVir = new AccumulatorAverageCovariance(blockSize);
        if (meterPrim != null) {
            DataPumpListener accumulatorPumpVir = new DataPumpListener(meterVir, accumulatorVir, interval);
            sim.integrator.getEventManager().addListener(accumulatorPumpVir);
        }

        //3 Centroid Virial
        AccumulatorAverageCovariance accumulatorCentVir = new AccumulatorAverageCovariance(blockSize);
        if (meterCentVir != null) {
            DataPumpListener accumulatorPumpCentVir = new DataPumpListener(meterCentVir, accumulatorCentVir, interval);
            sim.integrator.getEventManager().addListener(accumulatorPumpCentVir);
        }

         //4 HMAc (CLassical EC)
        AccumulatorAverageCovariance accumulatorHMAc = new AccumulatorAverageCovariance(blockSize);
        if (meterHMAc != null) {
            DataPumpListener accumulatorPumpHMAc = new DataPumpListener(meterHMAc, accumulatorHMAc, interval);
            sim.integrator.getEventManager().addListener(accumulatorPumpHMAc);
        }

        //5 HMAq (Quantum EC)
        AccumulatorAverageCovariance accumulatorNMSimple = new AccumulatorAverageCovariance(blockSize);
        if (meterNMSimple != null) {
            DataPumpListener accumulatorPumpHMAsimple = new DataPumpListener(meterNMSimple, accumulatorNMSimple, interval);
            sim.integrator.getEventManager().addListener(accumulatorPumpHMAsimple);
        }

        AccumulatorAverageCovariance accumulatorNMEC = new AccumulatorAverageCovariance(blockSize);
        if (meterNMEC != null) {
            DataPumpListener accumulatorPumpHMA = new DataPumpListener(meterNMEC, accumulatorNMEC, interval);
            sim.integrator.getEventManager().addListener(accumulatorPumpHMA);
        }

        AccumulatorAverageCovariance accumulatorStageSimple = new AccumulatorAverageCovariance(blockSize);
        if (meterStageSimple != null) {
            DataPumpListener pumpStagesimple = new DataPumpListener(meterStageSimple, accumulatorStageSimple, interval);
            sim.integrator.getEventManager().addListener(pumpStagesimple);
        }

        AccumulatorAverageCovariance accumulatorStageEC = new AccumulatorAverageCovariance(blockSize);
        if (meterStageEC != null) {
            DataPumpListener pumpStageEC = new DataPumpListener(meterStageEC, accumulatorStageEC, interval);
            sim.integrator.getEventManager().addListener(pumpStageEC);
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

        double kB_beta2 = sim.betaN*sim.betaN*nBeads*nBeads;

        //1 Prim
        DataGroup dataPrim = (DataGroup) accumulatorPrim.getData();
        IData dataAvgPrim = dataPrim.getData(accumulatorPrim.AVERAGE.index);
        IData dataErrPrim = dataPrim.getData(accumulatorPrim.ERROR.index);
        IData dataCorPrim = dataPrim.getData(accumulatorPrim.BLOCK_CORRELATION.index);
        IData dataCovPrim = dataPrim.getData(accumulatorPrim.COVARIANCE.index);

        double avgEnPrim = dataAvgPrim.getValue(0);
        double errEnPrim = dataErrPrim.getValue(0);
        double corEnPrim = dataCorPrim.getValue(0);
        System.out.println("\n En_prim:         " + avgEnPrim + "   err: " + errEnPrim + " cor: " + corEnPrim);
//        double CvnPrim = kB_beta2*(dataAvgPrim.getValue(1) - avgEnPrim*avgEnPrim);
//        double varX0 = errEnPrim*errEnPrim;
//        double varX1 = dataErrPrim.getValue(1)*dataErrPrim.getValue(1);
//        double corX0X1 = dataCovPrim.getValue(1)/Math.sqrt(dataCovPrim.getValue(0))/Math.sqrt(dataCovPrim.getValue(3));
//        double errCvnPrim = Math.sqrt(kB_beta2*(varX1 + 4.0*avgEnPrim*avgEnPrim*varX0 - 4*avgEnPrim*dataErrPrim.getValue(0)*dataErrPrim.getValue(1)*corX0X1));


        //2 Vir
        DataGroup dataVir = (DataGroup) accumulatorVir.getData();
        IData dataAvgVir = dataVir.getData(accumulatorVir.AVERAGE.index);
        IData dataErrVir = dataVir.getData(accumulatorVir.ERROR.index);
        IData dataCorVir = dataVir.getData(accumulatorVir.BLOCK_CORRELATION.index);
        IData dataCovVir = dataVir.getData(accumulatorVir.COVARIANCE.index);

        double avgEnVir = dataAvgVir.getValue(0);
        double errEnVir = dataErrVir.getValue(0);
        double corEnVir = dataCorVir.getValue(0);
        System.out.println(" En_vir:          " + avgEnVir + "   err: " + errEnVir + " cor: " + corEnVir);
//        double CvnVir = kB_beta2*(dataAvgVir.getValue(1) - avgEnVir*avgEnVir);
//        varX0 = errEnVir*errEnVir;
//        varX1 = dataErrVir.getValue(1)*dataErrVir.getValue(1);
//        corX0X1 = dataCovVir.getValue(1)/Math.sqrt(dataCovVir.getValue(0))/Math.sqrt(dataCovVir.getValue(3));
//        double errCvnVir = Math.sqrt(kB_beta2*(varX1 + 4.0*avgEnVir*avgEnVir*varX0 - 4*avgEnVir*dataErrVir.getValue(0)*dataErrVir.getValue(1)*corX0X1));

        //3 Cent-Vir
        DataGroup dataCentVir = (DataGroup) accumulatorCentVir.getData();
        IData dataAvgCentVir = dataCentVir.getData(accumulatorCentVir.AVERAGE.index);
        IData dataErrCentVir = dataCentVir.getData(accumulatorCentVir.ERROR.index);
        IData dataCorCentVir = dataCentVir.getData(accumulatorCentVir.BLOCK_CORRELATION.index);
        IData dataCovCentVir = dataCentVir.getData(accumulatorCentVir.COVARIANCE.index);
        double avgEnCentVir = dataAvgCentVir.getValue(0);
        double errEnCentVir = dataErrCentVir.getValue(0);
        double corEnCentVir = dataCorCentVir.getValue(0);
        System.out.println(" En_cvir:         " + avgEnCentVir + "   err: " + errEnCentVir + " cor: " + corEnCentVir);
//        double CvnCentVir = kB_beta2*(dataAvgCentVir.getValue(1) - avgEnCentVir*avgEnCentVir);
//        varX0 = errEnCentVir*errEnCentVir;
//        varX1 = dataErrCentVir.getValue(1)*dataErrCentVir.getValue(1);
//        corX0X1 = dataCovCentVir.getValue(1)/Math.sqrt(dataCovCentVir.getValue(0))/Math.sqrt(dataCovCentVir.getValue(3));
//        double errCvnCentVir = Math.sqrt(kB_beta2*(varX1 + 4.0*avgEnCentVir*avgEnCentVir*varX0 - 4*avgEnCentVir*dataErrCentVir.getValue(0)*dataErrCentVir.getValue(1)*corX0X1));

        //4 HMAc
        DataGroup dataHMAc = (DataGroup) accumulatorHMAc.getData();
        IData dataAvgHMAc = dataHMAc.getData(accumulatorHMAc.AVERAGE.index);
        IData dataErrHMAc = dataHMAc.getData(accumulatorHMAc.ERROR.index);
        IData dataCorHMAc = dataHMAc.getData(accumulatorHMAc.BLOCK_CORRELATION.index);
        IData dataCovHMAc = dataHMAc.getData(accumulatorHMAc.COVARIANCE.index);

        double avgEnHMAc = dataAvgHMAc.getValue(0) ;
        double errEnHMAc = dataErrHMAc.getValue(0);
        double corEnHMAc = dataCorHMAc.getValue(0);
        System.out.println(" En_hmac:         " + avgEnHMAc + "   err: " + errEnHMAc + " cor: " + corEnHMAc);
//        double CvnHMAc = kB_beta2*(dataAvgHMAc.getValue(1) - avgEnHMAc*avgEnHMAc);
//        varX0 = errEnHMAc*errEnHMAc;
//        varX1 = dataErrHMAc.getValue(1)*dataErrHMAc.getValue(1);
//        corX0X1 = dataCovHMAc.getValue(1)/Math.sqrt(dataCovHMAc.getValue(0))/Math.sqrt(dataCovHMAc.getValue(3));
//        double errCvnHMAc = Math.sqrt(kB_beta2*(varX1 + 4.0*avgEnHMAc*avgEnHMAc*varX0 - 4*avgEnHMAc*dataErrHMAc.getValue(0)*dataErrHMAc.getValue(1)*corX0X1));
//        if (errEnHMAc < 1e-10){
//            errCvnHMAc = Math.sqrt(kB_beta2*(varX1 + 4.0*avgEnHMAc*avgEnHMAc*varX0));
//        }


        //5 HMA simple NM
        DataGroup dataNMsimple = (DataGroup) accumulatorNMSimple.getData();
        IData dataAvgNMsimple = dataNMsimple.getData(accumulatorNMSimple.AVERAGE.index);
        IData dataErrNMsimple = dataNMsimple.getData(accumulatorNMSimple.ERROR.index);
        IData dataCorNMsimple = dataNMsimple.getData(accumulatorNMSimple.BLOCK_CORRELATION.index);
        IData dataCovNMsimple = dataNMsimple.getData(accumulatorNMSimple.COVARIANCE.index);

        double avgEnNMSimple = dataAvgNMsimple.getValue(0);
        double errEnNMSimple = dataErrNMsimple.getValue(0);
        double corEnNMSimple = dataCorNMsimple.getValue(0);
        System.out.println(" En_nm_simple:    " + avgEnNMSimple + "   err: " + errEnNMSimple + " cor: " + corEnNMSimple);

//        double Cvn_nm_simple  = kB_beta2*(dataAvgNMsimple.getValue(1) - avgEnNMSimple*avgEnNMSimple);
//        varX0 = errEnNMSimple*errEnNMSimple;
//        varX1 = dataErrNMsimple.getValue(1)*dataErrNMsimple.getValue(1);
//        corX0X1 = dataCovNMsimple.getValue(1)/Math.sqrt(dataCovNMsimple.getValue(0))/Math.sqrt(dataCovNMsimple.getValue(3));
//        double errCvnNMsimple = Math.sqrt(kB_beta2*(varX1 + 4.0*avgEnNMSimple*avgEnNMSimple*varX0 - 4*avgEnNMSimple*dataErrNMsimple.getValue(0)*dataErrNMsimple.getValue(1)*corX0X1));

        //6 HMA EC NM
        DataGroup dataNMEC = (DataGroup) accumulatorNMEC.getData();
        IData dataAvgNMEC = dataNMEC.getData(accumulatorNMEC.AVERAGE.index);
        IData dataErrNMEC = dataNMEC.getData(accumulatorNMEC.ERROR.index);
        IData dataCorNMEC = dataNMEC.getData(accumulatorNMEC.BLOCK_CORRELATION.index);
        IData dataCovNMEC = dataNMEC.getData(accumulatorNMEC.COVARIANCE.index);

        double avgEnNMEC = dataAvgNMEC.getValue(0);
        double errEnNMEC = dataErrNMEC.getValue(0);
        double corEnNMEC = dataCorNMEC.getValue(0);
        System.out.println(" En_nm_ec:        " + avgEnNMEC + "   err: " + errEnNMEC + " cor: " + corEnNMEC);

//        double CvnNMEC  = kB_beta2*(dataAvgNMEC.getValue(1) - avgEnNMEC*avgEnNMEC);
//        varX0 = errEnNMEC*errEnNMEC;
//        varX1 = dataErrNMEC.getValue(1)*dataErrNMEC.getValue(1);
//        corX0X1 = dataCovNMEC.getValue(1)/Math.sqrt(dataCovNMEC.getValue(0))/Math.sqrt(dataCovNMEC.getValue(3));
//        double errCvnNMEC = Math.sqrt(kB_beta2*(varX1 + 4.0*avgEnNMEC*avgEnNMEC*varX0 - 4*avgEnNMEC*errEnNMEC*dataErrNMEC.getValue(1)*corX0X1));
//        if (errEnNMEC < 1e-10){
//            errCvnNMEC = Math.sqrt(kB_beta2*(varX1 + 4.0*avgEnNMEC*avgEnNMEC*varX0));
//        }

        // 7 HMA simple stage
        DataGroup dataStageSimple = (DataGroup) accumulatorStageSimple.getData();
        IData dataAvgStageSimple = dataStageSimple.getData(accumulatorStageSimple.AVERAGE.index);
        IData dataErrStageSimple = dataStageSimple.getData(accumulatorStageSimple.ERROR.index);
        IData dataCorStageSimple = dataStageSimple.getData(accumulatorStageSimple.BLOCK_CORRELATION.index);
        IData dataCovStageSimple = dataStageSimple.getData(accumulatorStageSimple.COVARIANCE.index);
        double avgEnStageSimple = dataAvgStageSimple.getValue(0);
        double errEnStageSimple = dataErrStageSimple.getValue(0);
        double corEnStageSimple = dataCorStageSimple.getValue(0);
        System.out.println(" En_stage_simple: " + avgEnStageSimple + "   err: " + errEnStageSimple + " cor: " + corEnStageSimple);

//        double Cvn_stage_simple  = kB_beta2*(dataAvgStageSimple.getValue(1) - avgEnStageSimple*avgEnStageSimple);
//        varX0 = errEnStageSimple*errEnStageSimple;
//        varX1 = dataErrStageSimple.getValue(1)*dataErrStageSimple.getValue(1);
//        corX0X1 = dataCovStageSimple.getValue(1)/Math.sqrt(dataCovStageSimple.getValue(0))/Math.sqrt(dataCovStageSimple.getValue(3));
//        double errCvnStageSimple = Math.sqrt(kB_beta2*(varX1 + 4.0*avgEnStageSimple*avgEnStageSimple*varX0 - 4*avgEnStageSimple*dataErrStageSimple.getValue(0)*dataErrStageSimple.getValue(1)*corX0X1));

        //8 HMA EC stage
        DataGroup dataStageEC = (DataGroup) accumulatorStageEC.getData();
        IData dataAvgStageEC = dataStageEC.getData(accumulatorStageEC.AVERAGE.index);
        IData dataErrStageEC = dataStageEC.getData(accumulatorStageEC.ERROR.index);
        IData dataCorStageEC = dataStageEC.getData(accumulatorStageEC.BLOCK_CORRELATION.index);
        IData dataCovStageEC = dataStageEC.getData(accumulatorStageEC.COVARIANCE.index);
        double avgEnStageEC = dataAvgStageEC.getValue(0);
        double errEnStageEC = dataErrStageEC.getValue(0);
        double corEnStageEC = dataCorStageEC.getValue(0);
        System.out.println(" En_stage_ec:     " + avgEnStageEC + "   err: " + errEnStageEC + " cor: " + corEnStageEC);

//        double CvnStageEC  = kB_beta2*(dataAvgStageEC.getValue(1) - avgEnStageEC*avgEnStageEC);
//        varX0 = errEnStageEC*errEnStageEC;
//        varX1 = dataErrStageEC.getValue(1)*dataErrStageEC.getValue(1);
//        corX0X1 = dataCovStageEC.getValue(1)/Math.sqrt(dataCovStageEC.getValue(0))/Math.sqrt(dataCovStageEC.getValue(3));
//        double errCvnStageEC = Math.sqrt(kB_beta2*(varX1 + 4.0*avgEnStageEC*avgEnStageEC*varX0 - 4*avgEnStageEC*errEnStageEC*dataErrStageEC.getValue(1)*corX0X1));
//        if (errEnStageEC < 1e-10){
//            errCvnStageEC = Math.sqrt(kB_beta2*(varX1 + 4.0*avgEnStageEC*avgEnStageEC*varX0));
//        }

//        System.out.println("\n Cvn_prim: " + CvnPrim + " err: " + errCvnPrim);
//        System.out.println(" Cvn_vir: " + CvnVir + " err: " + errCvnVir);
//        System.out.println(" Cvn_cvir: " + CvnCentVir + " err: " + errCvnCentVir);
//        System.out.println(" Cvn_hmac: " + CvnHMAc + " err: " + errCvnHMAc);
//        System.out.println(" Cvn_nm_simple: " + Cvn_nm_simple + " err: " + errCvnNMsimple);
//        System.out.println(" Cvn_nm_ec: " + CvnNMEC + " err: " + errCvnNMEC);
//        System.out.println(" Cvn_stage_simple: " + Cvn_stage_simple + " err: " + errCvnStageSimple);
//        System.out.println(" Cvn_stage_ec: " + CvnStageEC + " err: " + errCvnStageEC);




//        System.out.println(" \n MSD_sim: " + avgMSD + " +/- " + errMSD + " cor: " + corMSD);
//        System.out.println(" MSDc: " + temperature/ mass/omega2);
//        System.out.println(" MSDq: " + hbar/mass/omega*(0.5+1.0/(Math.exp(hbar*omega/temperature)-1.0)));


        long endTime = System.currentTimeMillis();
        System.out.println("\n time (min): " + (endTime - startTime)/60.0/1000.0);
    }

    public enum MoveChoice {Real, NM, NMEC, Stage, StageEC};

    public static class OctaneParams extends ParameterBase {
        public double temperature = 1;
        public double hbar = 1;
        public double k2 = 1;
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