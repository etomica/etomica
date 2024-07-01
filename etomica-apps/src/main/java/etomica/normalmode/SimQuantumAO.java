/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.normalmode;

import etomica.action.activity.ActivityIntegrate;
import etomica.atom.AtomType;
import etomica.atom.IAtom;
import etomica.box.Box;
import etomica.config.ConformationLinear;
import etomica.data.*;
import etomica.data.types.DataGroup;
import etomica.graphics.ColorScheme;
import etomica.graphics.DisplayTextBox;
import etomica.graphics.SimulationGraphic;
import etomica.integrator.Integrator;
import etomica.integrator.IntegratorMC;
import etomica.integrator.mcmove.MCMoveAtom;
import etomica.integrator.mcmove.MCMoveBox;
import etomica.integrator.mcmove.MCMoveMolecule;
import etomica.integrator.mcmove.MCMoveMoleculeRotate;
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

public class SimQuantumAO extends Simulation {
    public Box box;
    public IntegratorMC integrator;
    public MCMoveBox move;
    public MCMoveHOReal2 moveStageSimple, moveStageEC;
    public MCMoveMolecule translateMove;
    public MCMoveMoleculeRotate rotateMove;
    public PotentialComputeField pcP1, pcP1EnTIA;
    public PotentialMasterBonding pmBonding;
    public PotentialCompute pmAgg;
//    public IPotential1 p1ah;
    public P1Anharmonic234 p1ah;
    public P1AnharmonicTIA p1ahUeff;
    public double betaN;
    public int nBeads;
    public double k2_kin;
    public SimQuantumAO(Space space, MoveChoice coordType, double mass, int nBeads, double temperature, double omega, double k3, double k4, boolean isTIA, double hbar) {
        super(space);
        SpeciesGeneral species = new SpeciesBuilder(space)
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
        double beta = 1.0/temperature;
        betaN = beta/nBeads;
        double omegaN = Math.sqrt(nBeads)/(hbar*beta);
        k2_kin = nBeads == 1 ? 0 : mass*omegaN*omegaN;
        double omega2 = omega*omega;

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
            p1ahUeff = new P1AnharmonicTIA(space, omega2, k4, nBeads, mass*omegaN*omegaN, facUeff);
            pcP1.setFieldPotential(species.getLeafType(), p1ahUeff);
        } else {
            p1ah = new P1Anharmonic234(space, omega2/nBeads, k3/nBeads, k4/nBeads);
            pcP1.setFieldPotential(species.getLeafType(), p1ah);
        }

        PotentialComputeAggregate.localStorageDefault = false;
        pmAgg = new PotentialComputeAggregate(pmBonding, pcP1);

        double facEn = 3.0;
        P1AnharmonicTIA p1ahEn = new P1AnharmonicTIA(space, mass*omega2, k4, nBeads, mass*omegaN*omegaN, facEn);
        pcP1EnTIA = new PotentialComputeField(getSpeciesManager(), box);
        pcP1EnTIA.setFieldPotential(species.getLeafType(), p1ahEn);

        integrator = new IntegratorMC(pmAgg, random, temperature, box);

        if (coordType == MoveChoice.Real) {
            move = new MCMoveAtom(random, pmAgg, box);
            integrator.getMoveManager().addMCMove(move);
            integrator.getMoveManager().setFrequency(move, nBeads);
        } else if (coordType == MoveChoice.NM) {
            move = new MCMoveHO(space, pmAgg, random, temperature, 0, box, hbar);
            integrator.getMoveManager().addMCMove(move);
        } else if (coordType == MoveChoice.NMEC) {
            move = new MCMoveHO(space, pmAgg, random, temperature, omega2, box, hbar);
            integrator.getMoveManager().addMCMove(move);
        } else if (coordType == MoveChoice.Stage) {
            move = new MCMoveHOReal2(space, pmAgg, random, temperature, 0, box, hbar);
            integrator.getMoveManager().addMCMove(move);
        } else if (coordType == MoveChoice.StageFast) {
            move = new MCMoveHOReal3(space, pmAgg, random, temperature, 0, box, hbar);
            integrator.getMoveManager().addMCMove(move);
        } else if (coordType == MoveChoice.StageEC) {
            move = new MCMoveHOReal2(space, pmAgg, random, temperature, omega2, box, hbar);
            integrator.getMoveManager().addMCMove(move);
        } else if (coordType == MoveChoice.StageFastEC) {
            move = new MCMoveHOReal3(space, pmAgg, random, temperature, omega2, box, hbar);
            integrator.getMoveManager().addMCMove(move);
        }
        else {
            throw new RuntimeException("Unknown move choice "+coordType);
        }

        if (coordType == MoveChoice.Real || coordType == MoveChoice.NM || coordType == MoveChoice.Stage) {
            translateMove = new MCMoveMolecule(random, pcP1, box);
            integrator.getMoveManager().addMCMove(translateMove);
            if (space.D() == 3 && coordType == MoveChoice.Real) {
                rotateMove = new MCMoveMoleculeRotate(random, pcP1, box);
                integrator.getMoveManager().addMCMove(rotateMove);
            }
        }

        moveStageEC = new MCMoveHOReal2(space, pmAgg, random, temperature, omega2, box, hbar);
        moveStageSimple = new MCMoveHOReal2(space, pmAgg, random, temperature, 0, box, hbar);
    }

    public Integrator getIntegrator() {
        return integrator;
    }

    public static void main(String[] args) {
        final long startTime = System.currentTimeMillis();
        OctaneParams params = new OctaneParams();
        if (args.length > 0) {
            ParseArgs.doParseArgs(params, args);
        }
        else {
            // custom parameters
            params.steps = 1000000;
            params.hbar = 1;
            params.temperature = 0.5;
            params.omega = 1;
            params.k3 = 0.1;
            params.k4 = 0.01;
            params.nShifts = 0;
            params.onlyCentroid = true;

//            params.coordType = MoveChoice.Real;
//            params.coordType = MoveChoice.NM;
//            params.coordType = MoveChoice.NMEC;
//            params.coordType = MoveChoice.Stage;
            params.coordType = MoveChoice.StageEC;

            params.nBeads = 64;
        }

        int nShifts = params.nShifts;
        double temperature = params.temperature;
        double hbar = params.hbar;
        double mass = params.mass;
        double omega = params.omega;
        double k3 = params.k3;
        double k4 = params.k4;
        boolean isGraphic = params.isGraphic;
        long steps = params.steps;
        long stepsEq = steps/10;
        boolean isTIA = params.isTIA;
        MoveChoice coordType = params.coordType;
        boolean zerok0 = params.zerok0;
        boolean onlyCentroid = params.onlyCentroid;
        double x = 1/temperature*hbar*omega;
        int nBeads = params.nBeads;
        if (nBeads == -1){
            nBeads = (int) (20*x); //20*x and 30*x are good for HO and AO, resp.
        }

        double omegaN = Math.sqrt(nBeads)*temperature/hbar;
        double omega2 = omega*omega;
        if (isTIA){
            omega2 = omega2*(1.0 + omega2/12.0/(nBeads*omegaN*omegaN));
        }
//        double actualOmega2 = omega2;
//        if (zerok0) omega2 = 0;

        final SimQuantumAO sim = new SimQuantumAO(Space1D.getInstance(), coordType, mass, nBeads, temperature, omega, k3, k4, isTIA, hbar);
        sim.integrator.reset();
        System.out.println(" PIMC-" + coordType);
        System.out.println(" mass: " + mass);
        System.out.println(" T: " + temperature);
        System.out.println(" hbar: " + hbar);
        System.out.println(" w: " + Math.sqrt(omega2));
        System.out.println(" wn: " + omegaN  + " , w/sqrt(n): " + Math.sqrt(omega2)/Math.sqrt(nBeads));
        System.out.println(" x = beta*hbar*w = " + hbar*omega/temperature);
        System.out.println(" nBeads: " + nBeads);
        System.out.println(" nShifts: "+ nShifts);
        System.out.println(" steps: " +  steps + " stepsEq: " + stepsEq);
        System.out.println(" omega: " + omega);
        System.out.println(" k3: " + k3);
        System.out.println(" k4: " + k4);
        System.out.println(" isTIA: " + isTIA);
        System.out.println(" onlyCentroid: " + onlyCentroid);

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
        double tau = 1/temperature*hbar*omega/nBeads;
        double D = 2 + omega2 / (nBeads*omegaN*omegaN);
        double alpha_2 = Math.log(D/2 + Math.sqrt(D*D/4 - 1));

        double EnQ_2 = temperature*(nBeads*Math.tanh(alpha_2/2.0)/Math.tanh(nBeads*alpha_2/2.0));
        System.out.println(" En_ho_q_2: " + EnQ_2);



        System.out.println(" E_ho_q: " + EnQinf);
        System.out.println(" Cvn_ho_c: " + CvnC);
        System.out.println(" Cvn_ho_q: " + CvnQ);
        System.out.println(" Cv_ho_c: " + CvnQinf + "\n");

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
//            meterHMAc = null;
//            meterHMA = new MeterPIHMATIA(sim.pmBonding, sim.pcP1EnTIA, sim.pcP1, sim.betaN, nBeads, omega2, sim.box);
//            meterHMAcent = null;
        } else {
            meterCentVir = new MeterPICentVir(sim.pcP1, temperature, nBeads, sim.box);
            if (!onlyCentroid) {
                meterPrim = new MeterPIPrim(sim.pmBonding, sim.pcP1, nBeads, temperature, sim.box);
                meterVir = new MeterPIVir(sim.pcP1, temperature, sim.box);
                meterHMAc = new MeterPIHMAc(sim.pcP1, temperature, nBeads, sim.box);
                meterNMSimple = new MeterPIHMA(sim.pmBonding, sim.pcP1, sim.betaN, nBeads, 0, sim.box, hbar);
                meterNMEC = new MeterPIHMA(sim.pmBonding, sim.pcP1, sim.betaN, nBeads, omega2, sim.box, hbar);
                meterStageSimple = new MeterPIHMAReal2(sim.pmBonding, sim.pcP1, nBeads, temperature, sim.moveStageSimple);
                meterStageSimple.setNumShifts(nShifts);
                meterStageEC = new MeterPIHMAReal2(sim.pmBonding, sim.pcP1, nBeads, temperature, sim.moveStageEC);
                meterStageEC.setNumShifts(nShifts);
            }
        }

//        MeterPIHMAvir meterHMAvir = new MeterPIHMAvir(sim.pmBonding, sim.pcP1, sim.betaN, nBeads, omega2, sim.box, hbar);//Bad!!

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
            DisplayTextBox timer = new DisplayTextBox();
            DataSourceCountSteps counter = new DataSourceCountSteps(sim.integrator);
            DataPumpListener counterPump = new DataPumpListener(counter, timer, 100);
            sim.integrator.getEventManager().addListener(counterPump);
            simGraphic.getPanel().controlPanel.add(timer.graphic());
            simGraphic.makeAndDisplayFrame("PIMC - "+coordType);
            return;
        }

        System.out.flush();


        //Write En
        DataLogger dlPrim = new DataLogger();
        DataLogger dlV = new DataLogger();
        DataLogger dlCV = new DataLogger();
        DataLogger dlHMAc = new DataLogger();
        DataLogger dlNMSimple = new DataLogger();
        DataLogger dlNMEC = new DataLogger();
        DataLogger dlStageSimple = new DataLogger();
        DataLogger dlStageEC = new DataLogger();
        boolean writeEn = false;
        if (writeEn) {
            int intervalDL = 1;
            int stepsDL = 10000;
            //Prim
            DataPumpListener dlPumpPrim = new DataPumpListener(meterPrim, dlPrim, intervalDL);
            sim.integrator.getEventManager().addListener(dlPumpPrim);
            dlPrim.setFileName("En_prim_T" + temperature + ".dat");
            dlPrim.setAppending(false);
            DataArrayWriter writerPrim = new DataArrayWriter();
            writerPrim.setIncludeHeader(false);
            dlPrim.setDataSink(writerPrim);

            //Vir
            DataPumpListener dlPumpV = new DataPumpListener(meterVir, dlV, intervalDL);
            sim.integrator.getEventManager().addListener(dlPumpV);
            dlV.setFileName("En_vir_T"+temperature+".dat");
            dlV.setAppending(false);
            DataArrayWriter writerV = new DataArrayWriter();
            writerV.setIncludeHeader(false);
            dlV.setDataSink(writerV);

            //CVir
            DataPumpListener dlPumpCV = new DataPumpListener(meterCentVir, dlCV, intervalDL);
            sim.integrator.getEventManager().addListener(dlPumpCV);
            dlCV.setFileName("En_cvir_T" + temperature + ".dat");
            dlCV.setAppending(false);
            DataArrayWriter writerCV = new DataArrayWriter();
            writerCV.setIncludeHeader(false);
            dlCV.setDataSink(writerCV);

            //HMAc
            DataPumpListener dlPumpHMAc = new DataPumpListener(meterHMAc, dlHMAc, intervalDL);
            sim.integrator.getEventManager().addListener(dlPumpHMAc);
            dlHMAc.setFileName("En_hmac_T" + temperature + ".dat");
            dlHMAc.setAppending(false);
            DataArrayWriter writerHMAc = new DataArrayWriter();
            writerHMAc.setIncludeHeader(false);
            dlHMAc.setDataSink(writerHMAc);

            //Simple NM
            DataPumpListener dlPumpNMSimple = new DataPumpListener(meterNMSimple, dlNMSimple, intervalDL);
            sim.integrator.getEventManager().addListener(dlPumpNMSimple);
            dlNMSimple.setFileName("En_nm_simple_T" + temperature + ".dat");
            dlNMSimple.setAppending(false);
            DataArrayWriter writerNMSimple = new DataArrayWriter();
            writerNMSimple.setIncludeHeader(false);
            dlNMSimple.setDataSink(writerNMSimple);

            // EC NM
            DataPumpListener dlPumpNMEC = new DataPumpListener(meterNMEC, dlNMEC, intervalDL);
            sim.integrator.getEventManager().addListener(dlPumpNMEC);
            dlNMEC.setFileName("En_nm_ec_T" + temperature + ".dat");
            dlNMEC.setAppending(false);
            DataArrayWriter writerNMEC = new DataArrayWriter();
            writerNMEC.setIncludeHeader(false);
            dlNMEC.setDataSink(writerNMEC);

            //Simple Stage
            DataPumpListener dlPumpStageSimple = new DataPumpListener(meterStageSimple, dlStageSimple, intervalDL);
            sim.integrator.getEventManager().addListener(dlPumpStageSimple);
            dlStageSimple.setFileName("En_stage_simple_T" + temperature + ".dat");
            dlStageSimple.setAppending(false);
            DataArrayWriter writerStageSimple = new DataArrayWriter();
            writerStageSimple.setIncludeHeader(false);
            dlStageSimple.setDataSink(writerStageSimple);

            // EC Stage
            DataPumpListener dlPumpStageEC = new DataPumpListener(meterStageEC, dlStageEC, intervalDL);
            sim.integrator.getEventManager().addListener(dlPumpStageEC);
            dlStageEC.setFileName("En_stage_ec_T" + temperature + ".dat");
            dlStageEC.setAppending(false);
            DataArrayWriter writerStageEC = new DataArrayWriter();
            writerStageEC.setIncludeHeader(false);
            dlStageEC.setDataSink(writerStageEC);

            sim.getController().runActivityBlocking(new ActivityIntegrate(sim.integrator, stepsDL));

            dlPrim.cleanUp();
            dlV.cleanUp();
            dlCV.cleanUp();
            dlHMAc.cleanUp();
            dlNMSimple.cleanUp();
            dlNMEC.cleanUp();
            dlStageSimple.cleanUp();
            dlStageEC.cleanUp();
            System.exit(0);
        }


        // equilibration
        sim.getController().runActivityBlocking(new ActivityIntegrate(sim.integrator, stepsEq));
        System.out.println(" equilibration finished");

        int interval;
        if (coordType == MoveChoice.Real) {
            interval = nBeads;
        } else {
            interval = 1;
        }
        int blocks = 100;
        long blockSize = steps/(interval*blocks);

        if (blockSize == 0) blockSize = 1;
        System.out.println(" blocks: " + blocks + " blocksize: " + blockSize + " interval: " + interval);

//        AccumulatorAverageFixed accumulatorMSD = new AccumulatorAverageFixed(blockSize);
//        DataPumpListener accumulatorPumpMSD = new DataPumpListener(meterMSDHO, accumulatorMSD, interval);
//        sim.integrator.getEventManager().addListener(accumulatorPumpMSD);

        //1 Primitive
        AccumulatorAverageCovariance accumulatorPrim = new AccumulatorAverageCovariance(blockSize);
        if (meterPrim != null) {
//            DataPumpListener accumulatorPumpPrim = new DataPumpListener(meterPrim, accumulatorPrim, interval);
//            sim.integrator.getEventManager().addListener(accumulatorPumpPrim);
        }

        //2 Virial
        AccumulatorAverageCovariance accumulatorVir = new AccumulatorAverageCovariance(blockSize);
        if (meterVir != null) {
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
        sim.integrator.getMoveManager().setEquilibrating(false);
        sim.getController().runActivityBlocking(new ActivityIntegrate(sim.integrator, steps));



        //MSD
//        DataGroup dataMSD = (DataGroup)accumulatorMSD.getData();
//        IData dataMSDAvg = dataMSD.getData(accumulatorMSD.AVERAGE.index);
//        IData dataMSDErr = dataMSD.getData(accumulatorMSD.ERROR.index);
//        IData dataMSDCorrelation = dataMSD.getData(accumulatorMSD.BLOCK_CORRELATION.index);
//        double avgMSD = dataMSDAvg.getValue(0);
//        double errMSD = dataMSDErr.getValue(0);
//        double corMSD = dataMSDCorrelation.getValue(0);


        double kB_beta2 = sim.betaN*sim.betaN*nBeads*nBeads;
        double varX0, varX1, corX0X1;



        // Cent-Vir
        DataGroup dataCentVir = (DataGroup) accumulatorCentVir.getData();
        IData dataAvgCentVir = dataCentVir.getData(accumulatorCentVir.AVERAGE.index);
        IData dataErrCentVir = dataCentVir.getData(accumulatorCentVir.ERROR.index);
        IData dataCorCentVir = dataCentVir.getData(accumulatorCentVir.BLOCK_CORRELATION.index);
        IData dataCovCentVir = dataCentVir.getData(accumulatorCentVir.COVARIANCE.index);
        double avgEnCentVir = dataAvgCentVir.getValue(0);
        double errEnCentVir = dataErrCentVir.getValue(0);
        double corEnCentVir = dataCorCentVir.getValue(0);
        double CvnCentVir = kB_beta2*(dataAvgCentVir.getValue(1) - avgEnCentVir*avgEnCentVir);
        varX0 = errEnCentVir*errEnCentVir;
        varX1 = dataErrCentVir.getValue(1)*dataErrCentVir.getValue(1);
        corX0X1 = dataCovCentVir.getValue(1)/Math.sqrt(dataCovCentVir.getValue(0))/Math.sqrt(dataCovCentVir.getValue(3));
        double errCvnCentVir = kB_beta2*Math.sqrt(varX1 + 4.0*avgEnCentVir*avgEnCentVir*varX0 - 4*avgEnCentVir*dataErrCentVir.getValue(0)*dataErrCentVir.getValue(1)*corX0X1);


        double avgEnPrim=0, avgEnVir=0, avgEnHMAc=0, avgEnNMSimple=0, avgEnNMEC=0, avgEnStageSimple=0,avgEnStageEC=0;
        double errEnPrim=0, errEnVir=0, errEnHMAc=0, errEnNMSimple=0, errEnNMEC=0, errEnStageSimple=0, errEnStageEC=0;
        double corEnPrim=0, corEnVir=0, corEnHMAc=0, corEnNMSimple=0, corEnNMEC=0, corEnStageSimple=0, corEnStageEC=0;
        double CvnPrim=0, CvnVir=0, CvnHMAc=0, CvnNMSimple=0, CvnNMEC=0, CvnStageSimple=0, CvnStageEC=0;
        double errCvnPrim=0, errCvnVir=0, errCvnHMAc=0, errCvnNMsimple=0, errCvnNMEC=0, errCvnStageSimple=0, errCvnStageEC=0;

        if (!onlyCentroid) {
            // Prim
            DataGroup dataPrim = (DataGroup) accumulatorPrim.getData();
            IData dataAvgPrim = dataPrim.getData(accumulatorPrim.AVERAGE.index);
            IData dataErrPrim = dataPrim.getData(accumulatorPrim.ERROR.index);
            IData dataCorPrim = dataPrim.getData(accumulatorPrim.BLOCK_CORRELATION.index);
            IData dataCovPrim = dataPrim.getData(accumulatorPrim.COVARIANCE.index);
            avgEnPrim = dataAvgPrim.getValue(0);
            errEnPrim = dataErrPrim.getValue(0);
            corEnPrim = dataCorPrim.getValue(0);
            CvnPrim = kB_beta2*(dataAvgPrim.getValue(1) - avgEnPrim*avgEnPrim);
            varX0 = errEnPrim*errEnPrim;
            varX1 = dataErrPrim.getValue(1)*dataErrPrim.getValue(1);
            corX0X1 = dataCovPrim.getValue(1)/Math.sqrt(dataCovPrim.getValue(0))/Math.sqrt(dataCovPrim.getValue(3));
            errCvnPrim = kB_beta2*Math.sqrt(varX1 + 4.0*avgEnPrim*avgEnPrim*varX0 - 4*avgEnPrim*dataErrPrim.getValue(0)*dataErrPrim.getValue(1)*corX0X1);

            // Vir
            DataGroup dataVir = (DataGroup) accumulatorVir.getData();
            IData dataAvgVir = dataVir.getData(accumulatorVir.AVERAGE.index);
            IData dataErrVir = dataVir.getData(accumulatorVir.ERROR.index);
            IData dataCorVir = dataVir.getData(accumulatorVir.BLOCK_CORRELATION.index);
            IData dataCovVir = dataVir.getData(accumulatorVir.COVARIANCE.index);
            avgEnVir = dataAvgVir.getValue(0);
            errEnVir = dataErrVir.getValue(0);
            corEnVir = dataCorVir.getValue(0);
            CvnVir = kB_beta2*(dataAvgVir.getValue(1) - avgEnVir*avgEnVir);
            varX0 = errEnVir*errEnVir;
            varX1 = dataErrVir.getValue(1)*dataErrVir.getValue(1);
            corX0X1 = dataCovVir.getValue(1)/Math.sqrt(dataCovVir.getValue(0))/Math.sqrt(dataCovVir.getValue(3));
            errCvnVir = kB_beta2*Math.sqrt(varX1 + 4.0*avgEnVir*avgEnVir*varX0 - 4*avgEnVir*dataErrVir.getValue(0)*dataErrVir.getValue(1)*corX0X1);

            // HMAc
            DataGroup dataHMAc = (DataGroup) accumulatorHMAc.getData();
            IData dataAvgHMAc = dataHMAc.getData(accumulatorHMAc.AVERAGE.index);
            IData dataErrHMAc = dataHMAc.getData(accumulatorHMAc.ERROR.index);
            IData dataCorHMAc = dataHMAc.getData(accumulatorHMAc.BLOCK_CORRELATION.index);
            IData dataCovHMAc = dataHMAc.getData(accumulatorHMAc.COVARIANCE.index);
            avgEnHMAc = dataAvgHMAc.getValue(0) ;
            errEnHMAc = dataErrHMAc.getValue(0);
            corEnHMAc = dataCorHMAc.getValue(0);
            CvnHMAc = kB_beta2*(dataAvgHMAc.getValue(1) - avgEnHMAc*avgEnHMAc);
            varX0 = errEnHMAc*errEnHMAc;
            varX1 = dataErrHMAc.getValue(1)*dataErrHMAc.getValue(1);
            corX0X1 = dataCovHMAc.getValue(1)/Math.sqrt(dataCovHMAc.getValue(0))/Math.sqrt(dataCovHMAc.getValue(3));
            errCvnHMAc = kB_beta2*Math.sqrt(varX1 + 4.0*avgEnHMAc*avgEnHMAc*varX0 - 4*avgEnHMAc*dataErrHMAc.getValue(0)*dataErrHMAc.getValue(1)*corX0X1);
            if (errEnHMAc < 1e-10){
                errCvnHMAc = kB_beta2*Math.sqrt(varX1 + 4.0*avgEnHMAc*avgEnHMAc*varX0 );
            }


            // HMA simple NM
            DataGroup dataNMsimple = (DataGroup) accumulatorNMSimple.getData();
            IData dataAvgNMsimple = dataNMsimple.getData(accumulatorNMSimple.AVERAGE.index);
            IData dataErrNMsimple = dataNMsimple.getData(accumulatorNMSimple.ERROR.index);
            IData dataCorNMsimple = dataNMsimple.getData(accumulatorNMSimple.BLOCK_CORRELATION.index);
            IData dataCovNMsimple = dataNMsimple.getData(accumulatorNMSimple.COVARIANCE.index);
            avgEnNMSimple = dataAvgNMsimple.getValue(0);
            errEnNMSimple = dataErrNMsimple.getValue(0);
            corEnNMSimple = dataCorNMsimple.getValue(0);
            CvnNMSimple  = kB_beta2*(dataAvgNMsimple.getValue(1) - avgEnNMSimple*avgEnNMSimple);
            varX0 = errEnNMSimple*errEnNMSimple;
            varX1 = dataErrNMsimple.getValue(1)*dataErrNMsimple.getValue(1);
            corX0X1 = dataCovNMsimple.getValue(1)/Math.sqrt(dataCovNMsimple.getValue(0))/Math.sqrt(dataCovNMsimple.getValue(3));
            errCvnNMsimple = kB_beta2*Math.sqrt(varX1 + 4.0*avgEnNMSimple*avgEnNMSimple*varX0 - 4*avgEnNMSimple*dataErrNMsimple.getValue(0)*dataErrNMsimple.getValue(1)*corX0X1);

            // HMA EC NM
            DataGroup dataNMEC = (DataGroup) accumulatorNMEC.getData();
            IData dataAvgNMEC = dataNMEC.getData(accumulatorNMEC.AVERAGE.index);
            IData dataErrNMEC = dataNMEC.getData(accumulatorNMEC.ERROR.index);
            IData dataCorNMEC = dataNMEC.getData(accumulatorNMEC.BLOCK_CORRELATION.index);
            IData dataCovNMEC = dataNMEC.getData(accumulatorNMEC.COVARIANCE.index);
            avgEnNMEC = dataAvgNMEC.getValue(0);
            errEnNMEC = dataErrNMEC.getValue(0);
            corEnNMEC = dataCorNMEC.getValue(0);
            CvnNMEC  = kB_beta2*(dataAvgNMEC.getValue(1) - avgEnNMEC*avgEnNMEC);
            varX0 = errEnNMEC*errEnNMEC;
            varX1 = dataErrNMEC.getValue(1)*dataErrNMEC.getValue(1);
            corX0X1 = dataCovNMEC.getValue(1)/Math.sqrt(dataCovNMEC.getValue(0))/Math.sqrt(dataCovNMEC.getValue(3));
            errCvnNMEC = kB_beta2*Math.sqrt(varX1 + 4.0*avgEnNMEC*avgEnNMEC*varX0 - 4*avgEnNMEC*errEnNMEC*dataErrNMEC.getValue(1)*corX0X1);
            if (errEnNMEC < 1e-10){
                errCvnNMEC = kB_beta2*Math.sqrt(varX1 + 4.0*avgEnNMEC*avgEnNMEC*varX0);
            }

            // HMA simple stage
            DataGroup dataStageSimple = (DataGroup) accumulatorStageSimple.getData();
            IData dataAvgStageSimple = dataStageSimple.getData(accumulatorStageSimple.AVERAGE.index);
            IData dataErrStageSimple = dataStageSimple.getData(accumulatorStageSimple.ERROR.index);
            IData dataCorStageSimple = dataStageSimple.getData(accumulatorStageSimple.BLOCK_CORRELATION.index);
            IData dataCovStageSimple = dataStageSimple.getData(accumulatorStageSimple.COVARIANCE.index);
            avgEnStageSimple = dataAvgStageSimple.getValue(0);
            errEnStageSimple = dataErrStageSimple.getValue(0);
            corEnStageSimple = dataCorStageSimple.getValue(0);
            CvnStageSimple  = kB_beta2*(dataAvgStageSimple.getValue(1) - avgEnStageSimple*avgEnStageSimple);
            varX0 = errEnStageSimple*errEnStageSimple;
            varX1 = dataErrStageSimple.getValue(1)*dataErrStageSimple.getValue(1);
            corX0X1 = dataCovStageSimple.getValue(1)/Math.sqrt(dataCovStageSimple.getValue(0))/Math.sqrt(dataCovStageSimple.getValue(3));
            errCvnStageSimple = kB_beta2*Math.sqrt(varX1 + 4.0*avgEnStageSimple*avgEnStageSimple*varX0 - 4*avgEnStageSimple*dataErrStageSimple.getValue(0)*dataErrStageSimple.getValue(1)*corX0X1);

            // HMA EC stage
            DataGroup dataStageEC = (DataGroup) accumulatorStageEC.getData();
            IData dataAvgStageEC = dataStageEC.getData(accumulatorStageEC.AVERAGE.index);
            IData dataErrStageEC = dataStageEC.getData(accumulatorStageEC.ERROR.index);
            IData dataCorStageEC = dataStageEC.getData(accumulatorStageEC.BLOCK_CORRELATION.index);
            IData dataCovStageEC = dataStageEC.getData(accumulatorStageEC.COVARIANCE.index);
            avgEnStageEC = dataAvgStageEC.getValue(0);
            errEnStageEC = dataErrStageEC.getValue(0);
            corEnStageEC = dataCorStageEC.getValue(0);
            CvnStageEC  = kB_beta2*(dataAvgStageEC.getValue(1) - avgEnStageEC*avgEnStageEC);
            varX0 = errEnStageEC*errEnStageEC;
            varX1 = dataErrStageEC.getValue(1)*dataErrStageEC.getValue(1);
            corX0X1 = dataCovStageEC.getValue(1)/Math.sqrt(dataCovStageEC.getValue(0))/Math.sqrt(dataCovStageEC.getValue(3));
            errCvnStageEC = kB_beta2*Math.sqrt(varX1 + 4.0*avgEnStageEC*avgEnStageEC*varX0 - 4*avgEnStageEC*errEnStageEC*dataErrStageEC.getValue(1)*corX0X1);
            if (errEnStageEC < 1e-10){
                errCvnStageEC = kB_beta2*Math.sqrt(varX1 + 4.0*avgEnStageEC*avgEnStageEC*varX0);
            }
        }

        if (onlyCentroid) {
            System.out.println();
            System.out.println(" En_cvir:          " + avgEnCentVir + "   err: " + errEnCentVir + " cor: " + corEnCentVir);
//            System.out.println(" Cvn_prim:         " + CvnPrim +          "   err: " + errCvnPrim);
            System.out.println(" Cvn_cvir:         " + CvnCentVir +       "   err: " + errCvnCentVir);
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
            System.out.println(" Cvn_prim:         " + CvnPrim +          "   err: " + errCvnPrim);
            System.out.println(" Cvn_vir:          " + CvnVir +           "   err: " + errCvnVir);
//            System.out.println(" Cvn_cvir:         " + CvnCentVir +       "   err: " + errCvnCentVir);
            System.out.println(" Cvn_hmac:         " + CvnHMAc +          "   err: " + errCvnHMAc);
            System.out.println(" Cvn_nm_simple:    " + CvnNMSimple +    "   err: " + errCvnNMsimple);
            System.out.println(" Cvn_nm_ec:        " + CvnNMEC +          "   err: " + errCvnNMEC);
            System.out.println(" Cvn_stage_simple: " + CvnStageSimple + "   err: " + errCvnStageSimple);
            System.out.println(" Cvn_stage_ec:     " + CvnStageEC +       "   err: " + errCvnStageEC);
        }




        //Acceptance ratio
        System.out.println("\n acceptance %: " + 100*sim.move.getTracker().acceptanceProbability());
        if (sim.translateMove!= null) {
            System.out.println(" translate step size: " + sim.translateMove.getStepSize());
            System.out.println(" acceptance % (translate):" + 100*sim.translateMove.getTracker().acceptanceProbability());
        }

        long endTime = System.currentTimeMillis();
        System.out.println("\n time: (min) " + (endTime - startTime)/60.0/1000.0);
    }

    public enum MoveChoice {Real, NM, NMEC, Stage, StageEC, StageFast, StageFastEC};

    public static class OctaneParams extends ParameterBase {
        public double temperature = 1.0;
        public double hbar = 0.1;
        public double mass = 1;
        public double omega = 1;
        public double k3 = 0.1;
        public double k4 = 0.01;
        public long steps = 1_000_000;
        public boolean isGraphic = false;
        public boolean isTIA = false;
        public boolean zerok0 = false;
        public MoveChoice coordType = MoveChoice.Real;
        public boolean onlyCentroid = !false;
        public int nBeads = -1;
        public int nShifts = 0;
    }
}