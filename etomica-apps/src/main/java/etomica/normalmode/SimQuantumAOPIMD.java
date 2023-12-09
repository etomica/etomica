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
import etomica.data.AccumulatorAverageCovariance;
import etomica.data.DataPumpListener;
import etomica.data.DataSourceCountTime;
import etomica.data.IData;
import etomica.data.types.DataGroup;
import etomica.graphics.ColorScheme;
import etomica.graphics.DisplayBox;
import etomica.graphics.DisplayTextBox;
import etomica.graphics.SimulationGraphic;
import etomica.integrator.Integrator;
import etomica.integrator.IntegratorLangevin;
import etomica.integrator.IntegratorMD;
import etomica.potential.P1Anharmonic;
import etomica.potential.P1AnharmonicTIA;
import etomica.potential.P2Harmonic;
import etomica.potential.PotentialMasterBonding;
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
    public P1Anharmonic p1ah;
    public P1AnharmonicTIA p1ahUeff;
    public double betaN;
    public double k2_kin;
    public int nBeads;
    public MCMoveHOReal2 moveStageSimple, moveStageEC;


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
            p1ah = new P1Anharmonic(space, k2, k4, nBeads);
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
            params.steps = 1000000;
            params.temperature = 1;
            params.k2 = 1;
            params.k4 = 24;
//            params.coordType = MoveChoice.Real;
//            params.coordType = MoveChoice.NM;
//            params.coordType = MoveChoice.NMEC;
//            params.coordType = MoveChoice.Stage;
            params.coordType = MoveChoice.StageEC;
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
        double w0 = Math.sqrt(k2/mass);
        double x = 1/temperature*hbar*w0;
        int nBeads = params.nBeads;
        if (nBeads == -1){
            nBeads = (int) (20*x);
        }

        double omegaN = Math.sqrt(nBeads)*temperature/hbar;
        double timeStep = params.timeStep;
        if (timeStep == -1) {
            double c = 0.1;
            if (params.coordType == MoveChoice.Real) {
                timeStep = c/omegaN/Math.sqrt(nBeads);//c/(20*w0)
            } else {
                timeStep = c/w0;
            }
        }

        double omega2 = k2/mass;
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
        System.out.println(" x = beta*hbar*w0: " + hbar*Math.sqrt(k2/mass)/temperature);
        System.out.println(" nBeads: " + nBeads);
        System.out.println(" nShifts: "+ nShifts);
        System.out.println(" timeStep: " + timeStep);
        System.out.println(" dt-non-real ~ 1/sqrt(omega2): " + 1/Math.sqrt(omega2));
        System.out.println(" dt-real ~ 1/omegaN: " + 1/omegaN/Math.sqrt(nBeads));

        System.out.println(" steps: " +  steps + " stepsEq: " + stepsEq);
        System.out.println(" k2: " + k2);
        System.out.println(" k4: " + k4);
        System.out.println(" isTIA: " + isTIA);
        System.out.println(" gammaLangevin: " + gammaLangevin);

        System.out.println("\n Quantum Harmonic Oscillator Theory");
        System.out.println(" ====================================");
        double omega = Math.sqrt(omega2);
        double alpha = 1 + 0.5*Math.pow(hbar*sim.betaN*omega,2)+0.5*hbar* sim.betaN*omega*Math.sqrt(4+Math.pow(hbar* sim.betaN*omega,2));
        double EnQ = sim.space.D()*(hbar*hbar*omega*omega)*sim.betaN*alpha/(alpha*alpha-1)*(Math.pow(alpha,nBeads)+1)/(Math.pow(alpha,nBeads)-1);

        double EnQinf = sim.space.D()*hbar*omega*(0.5 + 1/(Math.exp(nBeads*sim.betaN*hbar*omega)-1.0));
        double EnC = sim.space.D()*temperature;
        System.out.println(" EnC: " + EnC);
        System.out.println(" EnQ: " + EnQ);
        System.out.println(" EnQinf: " + EnQinf);

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
            meterCentVir = new MeterPICentVir(sim.pcP1, temperature, nBeads, sim.box);
            meterHMAc = new MeterPIHMAc(sim.pcP1, temperature, nBeads, sim.box);
            meterNMEC = new MeterPIHMA(sim.pmBonding, sim.pcP1, sim.betaN, nBeads, omega2, sim.box, hbar);
            meterNMSimple = new MeterPIHMA(sim.pmBonding, sim.pcP1, sim.betaN, nBeads, 0, sim.box, hbar);
            meterStageEC = new MeterPIHMAReal2(sim.pmBonding, sim.pcP1, nBeads, temperature, sim.moveStageEC);
            meterStageEC.setNumShifts(nShifts);
            meterStageSimple = new MeterPIHMAReal2(sim.pmBonding, sim.pcP1, nBeads, temperature, sim.moveStageSimple);
            meterStageSimple.setNumShifts(nShifts);
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
        System.out.println(" equilibration finished");

        int interval = 10;
        int blocks = 100;
        long blockSize = steps / (interval * blocks);
        if (blockSize == 0) blockSize = 1;
        System.out.println(" numBlocks: " + blocks + " blocksize: " + blockSize + " interval: " + interval);

//        AccumulatorAverageFixed accumulatorMSD = new AccumulatorAverageFixed(blockSize);
//        DataPumpListener accumulatorPumpMSD = new DataPumpListener(meterMSDHO, accumulatorMSD, interval);
//        sim.integrator.getEventManager().addListener(accumulatorPumpMSD);

        //1 Primitive
        AccumulatorAverageCovariance accumulatorPrim = new AccumulatorAverageCovariance(blockSize);
        if (meterPrim != null) {
            DataPumpListener accumulatorPumpPrim = new DataPumpListener(meterPrim, accumulatorPrim, interval);
            sim.integrator.getEventManager().addListener(accumulatorPumpPrim);
        }

        //2 Virial
        AccumulatorAverageCovariance accumulatorVir = new AccumulatorAverageCovariance(blockSize);
//        if (meterPrim != null) {
//            DataPumpListener accumulatorPumpVir = new DataPumpListener(meterVir, accumulatorVir, interval);
//            sim.integrator.getEventManager().addListener(accumulatorPumpVir);
//        }

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
        AccumulatorAverageCovariance accumulatorHMAsimple = new AccumulatorAverageCovariance(blockSize);
        if (meterNMSimple != null) {
            DataPumpListener accumulatorPumpHMAsimple = new DataPumpListener(meterNMSimple, accumulatorHMAsimple, interval);
            sim.integrator.getEventManager().addListener(accumulatorPumpHMAsimple);
        }

        AccumulatorAverageCovariance accumulatorHMA = new AccumulatorAverageCovariance(blockSize);
        if (meterNMEC != null) {
            DataPumpListener accumulatorPumpHMA = new DataPumpListener(meterNMEC, accumulatorHMA, interval);
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
        double kB_beta2 = sim.betaN*sim.betaN*nBeads*nBeads;




        //1 Prim
        if (meterPrim!=null) {
            DataGroup dataPrim = (DataGroup) accumulatorPrim.getData();
            IData dataAvgPrim = dataPrim.getData(accumulatorPrim.AVERAGE.index);
            IData dataErrPrim = dataPrim.getData(accumulatorPrim.ERROR.index);
            IData dataCorPrim = dataPrim.getData(accumulatorPrim.BLOCK_CORRELATION.index);
            IData dataCovPrim = dataPrim.getData(accumulatorPrim.COVARIANCE.index);

            double avgEnPrim = dataAvgPrim.getValue(0);
            double errEnPrim = dataErrPrim.getValue(0);
            double corEnPrim = dataCorPrim.getValue(0);
            System.out.println("\n En_prim:         " + avgEnPrim + "   err: " + errEnPrim + " cor: " + corEnPrim);
            double CvnPrim = kB_beta2*(dataAvgPrim.getValue(1) + dataCovPrim.getValue(0));
            System.out.println(" Cvn_prim: " + CvnPrim);
        }


        //2 Vir
        if (meterVir!=null) {
            DataGroup dataVir = (DataGroup) accumulatorVir.getData();
            IData dataAvgVir = dataVir.getData(accumulatorVir.AVERAGE.index);
            IData dataErrVir = dataVir.getData(accumulatorVir.ERROR.index);
            IData dataCorVir = dataVir.getData(accumulatorVir.BLOCK_CORRELATION.index);
            IData dataCovVir = dataVir.getData(accumulatorVir.COVARIANCE.index);

            double avgEnVir = dataAvgVir.getValue(0);
            double errEnVir = dataErrVir.getValue(0);
            double corEnVir = dataCorVir.getValue(0);
            System.out.println(" En_vir:          " + avgEnVir + "   err: " + errEnVir + " cor: " + corEnVir);
            double CvnVir = kB_beta2*(dataAvgVir.getValue(1) + dataCovVir.getValue(0));
            System.out.println(" Cvn_vir: " + CvnVir);
        }

        //3 Cent-Vir
        if (meterCentVir!=null) {
            DataGroup dataCentVir = (DataGroup) accumulatorCentVir.getData();
            IData dataAvgCentVir = dataCentVir.getData(accumulatorCentVir.AVERAGE.index);
            IData dataErrCentVir = dataCentVir.getData(accumulatorCentVir.ERROR.index);
            IData dataCorCentVir = dataCentVir.getData(accumulatorCentVir.BLOCK_CORRELATION.index);
            IData dataCovCentVir = dataCentVir.getData(accumulatorCentVir.COVARIANCE.index);

            double avgEnCentVir = dataAvgCentVir.getValue(0);
            double errEnCentVir = dataErrCentVir.getValue(0);
            double corEnCentVir = dataCorCentVir.getValue(0);
            System.out.println(" En_cvir:         " + avgEnCentVir + "   err: " + errEnCentVir + " cor: " + corEnCentVir);
            double CvnCentVir = kB_beta2*(dataAvgCentVir.getValue(1) + dataCovCentVir.getValue(0));
            System.out.println(" Cvn_cvir: " + CvnCentVir);
        }

        //4 HMAc
        if (meterHMAc!=null) {
            DataGroup dataHMAc = (DataGroup) accumulatorHMAc.getData();
            IData dataAvgHMAc = dataHMAc.getData(accumulatorHMAc.AVERAGE.index);
            IData dataErrHMAc = dataHMAc.getData(accumulatorHMAc.ERROR.index);
            IData dataCorHMAc = dataHMAc.getData(accumulatorHMAc.BLOCK_CORRELATION.index);
            IData dataCovHMAc = dataHMAc.getData(accumulatorHMAc.COVARIANCE.index);

            double avgEnHMAc = dataAvgHMAc.getValue(0) ;
            double errEnHMAc = dataErrHMAc.getValue(0);
            double corEnHMAc = dataCorHMAc.getValue(0);
            System.out.println(" En_hmac:         " + avgEnHMAc + "   err: " + errEnHMAc + " cor: " + corEnHMAc);
            double CvnHMAc = kB_beta2*(dataAvgHMAc.getValue(1) + dataCovHMAc.getValue(0));
            System.out.println(" Cvn_hmac: " + CvnHMAc);
        }


        //5 HMA simple NM
        if (meterNMSimple!=null) {
            DataGroup dataHMAsimple = (DataGroup) accumulatorHMAsimple.getData();
            IData dataAvgHMAsimple = dataHMAsimple.getData(accumulatorHMAsimple.AVERAGE.index);
            IData dataErrHMAsimple = dataHMAsimple.getData(accumulatorHMAsimple.ERROR.index);
            IData dataCorHMAsimple = dataHMAsimple.getData(accumulatorHMAsimple.BLOCK_CORRELATION.index);
            IData dataCovHMAsimple = dataHMAsimple.getData(accumulatorHMAsimple.COVARIANCE.index);

            double avgEnHMAsimple = dataAvgHMAsimple.getValue(0);
            double errEnHMAsimple = dataErrHMAsimple.getValue(0);
            double corEnHMAsimple = dataCorHMAsimple.getValue(0);
            System.out.println(" En_nm_simple:    " + avgEnHMAsimple + "   err: " + errEnHMAsimple + " cor: " + corEnHMAsimple);
            double Cvn_nm_simple  = kB_beta2*(dataAvgHMAsimple.getValue(1) + dataCovHMAsimple.getValue(0));
            System.out.println(" Cvn_nm_simple: " + Cvn_nm_simple);
        }

        //6 HMA EC NM
        if (meterNMEC!=null) {
            DataGroup dataHMA = (DataGroup) accumulatorHMA.getData();
            IData dataAvgHMA = dataHMA.getData(accumulatorHMA.AVERAGE.index);
            IData dataErrHMA = dataHMA.getData(accumulatorHMA.ERROR.index);
            IData dataCorHMA = dataHMA.getData(accumulatorHMA.BLOCK_CORRELATION.index);
            IData dataCovHMA = dataHMA.getData(accumulatorHMA.COVARIANCE.index);

            double avgEnHMA = dataAvgHMA.getValue(0);
            double errEnHMA = dataErrHMA.getValue(0);
            double corEnHMA = dataCorHMA.getValue(0);
            System.out.println(" En_nm_EC:        " + avgEnHMA + "   err: " + errEnHMA + " cor: " + corEnHMA);
            double CvnHMA  = kB_beta2*(dataAvgHMA.getValue(1) + dataCovHMA.getValue(0));
            System.out.println(" Cvn_hma: " + CvnHMA);
        }


        // 7 HMA simple stage
        if (meterStageSimple != null) {
            DataGroup dataStageECsimple = (DataGroup) accumulatorStageSimple.getData();
            IData dataAvgStageECsimple = dataStageECsimple.getData(accumulatorStageSimple.AVERAGE.index);
            IData dataErrStageECsimple = dataStageECsimple.getData(accumulatorStageSimple.ERROR.index);
            IData dataCorStageECsimple = dataStageECsimple.getData(accumulatorStageSimple.BLOCK_CORRELATION.index);
            IData dataCovStageECsimple = dataStageECsimple.getData(accumulatorStageSimple.COVARIANCE.index);
            double avgStageECsimple = dataAvgStageECsimple.getValue(0);
            double errStageECsimple = dataErrStageECsimple.getValue(0);
            double corStageECsimple = dataCorStageECsimple.getValue(0);
            System.out.println(" En_stage_simple: " + avgStageECsimple + "   err: " + errStageECsimple + " cor: " + corStageECsimple);
            double Cvn_stage_simple  = kB_beta2*(dataAvgStageECsimple.getValue(1) + dataCovStageECsimple.getValue(0));
            System.out.println(" Cvn_stage_simple: " + Cvn_stage_simple);
        }

        //8 HMA EC stage
        if (meterStageEC != null) {
            DataGroup dataStageEC = (DataGroup) accumulatorStageEC.getData();
            IData dataAvgStageEC = dataStageEC.getData(accumulatorStageEC.AVERAGE.index);
            IData dataErrStageEC = dataStageEC.getData(accumulatorStageEC.ERROR.index);
            IData dataCorStageEC = dataStageEC.getData(accumulatorStageEC.BLOCK_CORRELATION.index);
            IData dataCovStageEC = dataStageEC.getData(accumulatorStageEC.COVARIANCE.index);
            double avgStageEC = dataAvgStageEC.getValue(0);
            double errStageEC = dataErrStageEC.getValue(0);
            double corStageEC = dataCorStageEC.getValue(0);
            System.out.println(" En_stage_EC:     " + avgStageEC + "   err: " + errStageEC + " cor: " + corStageEC);
            double Cvn_stage_EC  = kB_beta2*(dataAvgStageEC.getValue(1) + dataCovStageEC.getValue(0));
            System.out.println(" Cvn_stage_EC: " + Cvn_stage_EC);
        }



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