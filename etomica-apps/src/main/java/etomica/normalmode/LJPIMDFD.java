/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.normalmode;

import etomica.action.BoxInflate;
import etomica.action.activity.ActivityIntegrate;
import etomica.atom.AtomType;
import etomica.atom.DiameterHashByType;
import etomica.atom.IAtom;
import etomica.box.Box;
import etomica.config.ConfigurationLattice;
import etomica.config.ConformationLinear;
import etomica.data.*;
import etomica.data.history.HistoryCollapsingAverage;
import etomica.data.meter.MeterPotentialEnergyFromIntegrator;
import etomica.data.types.DataFunction;
import etomica.data.types.DataGroup;
import etomica.graphics.*;
import etomica.integrator.IntegratorLangevin;
import etomica.integrator.IntegratorMD;
import etomica.lattice.LatticeCubicFcc;
import etomica.potential.*;
import etomica.potential.compute.PotentialComputeAggregate;
import etomica.potential.compute.PotentialComputePair;
import etomica.simulation.Simulation;
import etomica.space.Space;
import etomica.space3d.Space3D;
import etomica.species.SpeciesBuilder;
import etomica.species.SpeciesGeneral;
import etomica.util.ParameterBase;
import etomica.util.ParseArgs;

import java.awt.*;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

/**
 * PIMD using Langevin Thermostat using BAOAB algorithm.
 *
 * @author Andrew Schultz
 */

public class LJPIMDFD extends Simulation {

    public final PotentialComputePair potentialMaster;
    public final IntegratorMD integrator;
    public final PotentialMasterBonding pmBonding;
    public final Box box;
    public final PotentialComputeAggregate pmAgg;
    public final MCMoveHOReal2 moveStageSimple, moveStageEC;
    public double betaN;
    public int dim;

    public LJPIMDFD(Space space, MoveChoice coordType, double mass, double timeStep, double gammaLangevin, int nBeads, int numAtoms, double temperature, double density, double rc, double omega2, double hbar) {
        super(Space3D.getInstance());
        SpeciesGeneral species = new SpeciesBuilder(space)
                .setDynamic(true)
                .addCount(AtomType.simple("A", mass / nBeads), nBeads)
                .withConformation(new ConformationLinear(space, 0))
                .build();
        addSpecies(species);
        this.dim = space.D();
        box = new Box(space);
        addBox(box);
        NeighborListManagerPI neighborManager = new NeighborListManagerPI(getSpeciesManager(), box, 2, 1.4*rc, BondingInfo.noBonding());
        potentialMaster = new PotentialComputePair(getSpeciesManager(), box, neighborManager);
        pmBonding = new PotentialMasterBonding(getSpeciesManager(), box);
        double beta = 1 / temperature;
        betaN = beta/nBeads;

        double omegaN = nBeads/(hbar*beta);
        double k2_kin = nBeads == 1 ? 0 : mass*omegaN*omegaN/nBeads;
        P2Harmonic p2Bond = new P2Harmonic(k2_kin, 0);
        List<int[]> pairs = new ArrayList<>();
        for (int i = 0; i < nBeads; i++) {
            int[] p = new int[]{i, (i + 1) % nBeads};
            pairs.add(p);
        }
        pmBonding.setBondingPotentialPair(species, p2Bond, pairs);

        box.setNMolecules(species, numAtoms);
        BoxInflate inflater = new BoxInflate(box, space, density);
        inflater.actionPerformed();
        ConfigurationLattice config = new ConfigurationLattice(new LatticeCubicFcc(space), space);
        config.initializeCoordinates(box);

        double sigma = 1.0, epsilon = 1.0;
        P2LennardJones p2lj = new P2LennardJones(1.0, 1.0 / nBeads);
        P2SoftSphericalTruncatedForceShifted p2 = new P2SoftSphericalTruncatedForceShifted(p2lj, rc);
        AtomType atomType = species.getLeafType();
        potentialMaster.setPairPotential(atomType, atomType, p2);
        potentialMaster.doAllTruncationCorrection = false;

        PotentialComputeAggregate.localStorageDefault = true;
        pmAgg = new PotentialComputeAggregate(pmBonding, potentialMaster);

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

        moveStageSimple = new MCMoveHOReal2(space, pmAgg, random, temperature, 0, box, hbar);
        moveStageEC = new MCMoveHOReal2(space, pmAgg, random, temperature, omega2, box, hbar);
        integrator.setThermostatNoDrift(false);
        integrator.setIsothermal(true);

    }

    public static void main(String[] args) {
        long t1 = System.currentTimeMillis();
        SimParams params = new SimParams();
        if (args.length > 0) {
            ParseArgs.doParseArgs(params, args);
        } else {
            params.steps = 1000;
            params.hbar = 0.1;
            params.temperature = 1;
            params.numAtoms = 32;
            params.rc = 2.5;
            params.isGraphic = false;
            params.onlyCentroid = true;
//            params.coordType = MoveChoice.Real;
//            params.coordType = MoveChoice.NM;
//            params.coordType = MoveChoice.Stage;
//            params.coordType = MoveChoice.NMEC;
            params.coordType = MoveChoice.StageEC;
            params.nBeads = 16;
            params.timeStep = 0.01;
        }
        double facTimestep = params.facTimestep;
        Space space = Space.getInstance(params.D);
        int nShifts = params.nShifts;
        double mass = params.mass;
        double hbar = params.hbar;
        int numAtoms = params.numAtoms;
        int nBeads = params.nBeads;
        double temperature = params.temperature;
        double density = params.density;
        double rc = params.rc;
        double k2 = params.k2;
        double omega2 = params.k2 / mass;
        double gammaLangevin = params.gammaLangevin;
        double timeStep = params.timeStep;
        boolean isGraphic = params.isGraphic;
        boolean onlyCentroid = params.onlyCentroid;
        MoveChoice coordType = params.coordType;
        long steps = params.steps;
        long stepsEq = steps/10;

        double omega = Math.sqrt(omega2);
        double x = 1/temperature*hbar*omega;
        if (nBeads == -1){
            nBeads = (int) (20*x);
        }

        double omegaN = nBeads*temperature/hbar;
//        if (timeStep == -1) {
//            if (coordType == MoveChoice.Real) {
//                timeStep = facTimestep/omegaN/Math.sqrt(nBeads);// mi=m/n in real space, so m wn^2 = (m/n)*(n wn^2)==> dt~1/[wn sqrt(n)]
//            } else if (coordType == MoveChoice.NM) {
//                double s = omega2/omegaN/omegaN/nBeads;
//                timeStep = facTimestep*Math.sqrt(s)/omega; // which is 1/sqrt(n)
//            } else if (coordType == MoveChoice.Stage) {
//                double s = omega2/omegaN/omegaN * (1 + 1.0 / 12.0 * (nBeads * nBeads - 1.0) / nBeads);
//                timeStep = facTimestep*Math.sqrt(s)/omega; // for large n, timeStep ~ hbar/T
//            } else {
//                timeStep = facTimestep/omega;
//            }
//        }

        LJPIMDFD sim = new LJPIMDFD(space, coordType, mass, timeStep, gammaLangevin, nBeads, numAtoms, temperature, density, rc, omega2, hbar);
        sim.integrator.reset();

        System.out.println(" LJ PIMD-"+coordType);
        System.out.println(" mass: " + mass);
        System.out.println(" T: " + temperature);
        System.out.println(" hbar: " + hbar);
        System.out.println(" w: " + omega);
        System.out.println(" wn: " + omegaN  + " , w/sqrt(n): " + omega/Math.sqrt(nBeads));
        System.out.println(" x: beta*hbar*w = " + hbar*omega/temperature);
        System.out.println(" nBeads: " + nBeads);
        System.out.println(" nShifts: "+ nShifts);
        System.out.println(" steps: " +  steps + " stepsEq: " + stepsEq);
        System.out.println(" timestep: " + timeStep);
        System.out.println(" facTimestep: " + facTimestep);
        System.out.println(" k2: " + k2);
        System.out.println(" gammaLangevin: " + gammaLangevin);
        System.out.println(" N: " + numAtoms);
        System.out.println(" density: " + density);
        System.out.println(" rc: " + rc);

        DataSourceScalar meterKE = sim.integrator.getMeterKineticEnergy();

        MeterPIPrim meterPrim = null;
        MeterPIVir meterVir = null;
        MeterPICentVir meterCentVir = null;
        MeterPICentVirFD meterCentVirFD = null;
        MeterPIHMAc meterHMAc = null;
        MeterPIHMAcFD meterHMAcFD = null;
        MeterPIHMA meterNMEC = null;
        MeterPIHMA meterNMSimple = null;
        MeterPIHMAReal2 meterStageEC = null;
        MeterPIHMAReal2 meterStageSimple = null;

        meterCentVir = new MeterPICentVir(sim.potentialMaster, temperature, nBeads, sim.box);
        double dbeta = 0.01;
        meterCentVirFD = new MeterPICentVirFD(sim.potentialMaster, temperature, nBeads, sim.box, dbeta);
        meterPrim = new MeterPIPrim(sim.pmBonding, sim.potentialMaster, nBeads, temperature, sim.box);
        meterHMAc = new MeterPIHMAc(sim.potentialMaster, temperature, nBeads, sim.box);
        meterHMAcFD = new MeterPIHMAcFD(sim.potentialMaster, temperature, nBeads, sim.box, dbeta);

        if (!onlyCentroid) {
            meterPrim = new MeterPIPrim(sim.pmBonding, sim.potentialMaster, nBeads, temperature, sim.box);
            meterVir = new MeterPIVir(sim.potentialMaster, temperature, sim.box);
            meterCentVir = new MeterPICentVir(sim.potentialMaster, temperature, nBeads, sim.box);
            meterNMSimple = new MeterPIHMA(sim.pmBonding, sim.potentialMaster, sim.betaN, nBeads, 0, sim.box, hbar);
            meterNMEC = new MeterPIHMA(sim.pmBonding, sim.potentialMaster, sim.betaN, nBeads, omega2, sim.box, hbar);
            meterStageSimple = new MeterPIHMAReal2(sim.pmBonding, sim.potentialMaster, nBeads, temperature, sim.moveStageSimple, nShifts);
            meterStageEC = new MeterPIHMAReal2(sim.pmBonding, sim.potentialMaster, nBeads, temperature, sim.moveStageEC, nShifts);
        }


        int interval = 5;
        int blocks = 100;
        long blockSize = steps / (interval * blocks);
        if (blockSize == 0) blockSize = 1;
        System.out.println(" numBlocks: " + blocks + " blocksize: " + blockSize + " interval: " + interval);

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

            ((DiameterHashByType) ((DisplayBox) simGraphic.displayList().getFirst()).getDiameterHash()).setDiameter(sim.species().getAtomType(0), 0.2);

            DisplayTextBox timer = new DisplayTextBox();
            DataSourceCountTime counter = new DataSourceCountTime(sim.integrator);
            DataPumpListener counterPump = new DataPumpListener(counter, timer, 100);
            sim.integrator.getEventManager().addListener(counterPump);
            simGraphic.getPanel().controlPanel.add(timer.graphic());

            MeterMSDHO meterMSD = new MeterMSDHO(sim.box);
            AccumulatorHistory historyMSD = new AccumulatorHistory(new HistoryCollapsingAverage());
            historyMSD.setTimeDataSource(counter);
            DataPumpListener pumpMSD = new DataPumpListener(meterMSD, historyMSD, interval);
            sim.integrator.getEventManager().addListener(pumpMSD);
            DisplayPlotXChart plotMSD = new DisplayPlotXChart();
            plotMSD.setLabel("MSD");
            plotMSD.setDoLegend(false);
            historyMSD.addDataSink(plotMSD.makeSink("MSD history"));
            simGraphic.add(plotMSD);

            MeterPotentialEnergyFromIntegrator meterPE = new MeterPotentialEnergyFromIntegrator(sim.integrator);
            AccumulatorHistory historyPE = new AccumulatorHistory(new HistoryCollapsingAverage());
            historyPE.setTimeDataSource(counter);
            DataPumpListener pumpPE = new DataPumpListener(meterPE, historyPE, interval);
            sim.integrator.getEventManager().addListener(pumpPE);


            simGraphic.makeAndDisplayFrame("PIMD-"+coordType);

            return;
        }

        System.out.flush();

//        ConfigurationStorage configStorage = new ConfigurationStorage(sim.box, ConfigurationStorage.StorageType.MSD);
//        DataSourceBAC meterBAC = new DataSourceBAC(configStorage);
//        DataSourceRAC meterRAC = new DataSourceRAC(configStorage);
//        DataSourceMSDAC meterMSDAC = new DataSourceMSDAC(configStorage);

        //Write En
        DataLogger dlPrim = new DataLogger();
//        DataLogger dlV = new DataLogger();
        DataLogger dlCV = new DataLogger();
        DataLogger dlHMAc = new DataLogger();
        DataLogger dlHMA2 = new DataLogger();

        boolean writeEn = false;
        if (writeEn) {
            int intervalDL = 1;
            //Prim
            DataPumpListener dlPumpPrim = new DataPumpListener(meterPrim, dlPrim, intervalDL);
            sim.integrator.getEventManager().addListener(dlPumpPrim);
            dlPrim.setFileName("En_primT" + temperature + ".dat");
            dlPrim.setAppending(false);
            DataArrayWriter writerPrim = new DataArrayWriter();
            writerPrim.setIncludeHeader(false);
            dlPrim.setDataSink(writerPrim);
            //Vir
//            DataPumpListener dlPumpV = new DataPumpListener(meterVir, dlV, intervalDL);
//            sim.integrator.getEventManager().addListener(dlPumpV);
//            dlV.setFileName("En_virT"+temperature+".dat");
//            dlV.setAppending(false);
//            DataArrayWriter writerV = new DataArrayWriter();
//            writerV.setIncludeHeader(false);
//            dlV.setDataSink(writerV);
//
            //CVir
            DataPumpListener dlPumpCV = new DataPumpListener(meterCentVir, dlCV, intervalDL);
            sim.integrator.getEventManager().addListener(dlPumpCV);
            dlCV.setFileName("En_cvirT" + temperature + ".dat");
            dlCV.setAppending(false);
            DataArrayWriter writerCV = new DataArrayWriter();
            writerCV.setIncludeHeader(false);
            dlCV.setDataSink(writerCV);

            //HMAc
            DataPumpListener dlPumpHMAc = new DataPumpListener(meterHMAc, dlHMAc, intervalDL);
            sim.integrator.getEventManager().addListener(dlPumpHMAc);
            dlHMAc.setFileName("En_hmacT" + temperature + ".dat");
            dlHMAc.setAppending(false);
            DataArrayWriter writerHMAc = new DataArrayWriter();
            writerHMAc.setIncludeHeader(false);
            dlHMAc.setDataSink(writerHMAc);

            //HMA2
            DataPumpListener dlPumpHMA2 = new DataPumpListener(meterStageEC, dlHMA2, intervalDL);
            sim.integrator.getEventManager().addListener(dlPumpHMA2);
            dlHMA2.setFileName("En_hma2T" + temperature + ".dat");
            dlHMA2.setAppending(false);
            DataArrayWriter writerHMA2 = new DataArrayWriter();
            writerHMA2.setIncludeHeader(false);
            dlHMA2.setDataSink(writerHMA2);

            sim.getController().runActivityBlocking(new ActivityIntegrate(sim.integrator, 1000));
            dlPrim.cleanUp();
//            dlV.cleanUp();
            dlCV.cleanUp();
            dlHMAc.cleanUp();
            dlHMA2.cleanUp();
            System.exit(0);
        }

        sim.potentialMaster.computeAll(true);
        double uLat = sim.potentialMaster.getLastEnergy()/numAtoms;
        double volume = sim.box.getBoundary().volume();
        double pLat = -sim.potentialMaster.getLastVirial()/3.0/volume;
        System.out.println(" uLat: " + uLat);
        System.out.println(" pLat: " + pLat);


        // equilibration
        sim.getController().runActivityBlocking(new ActivityIntegrate(sim.integrator, stepsEq));
        System.out.println(" equilibration finished");

//        configStorage.setEnabled(true);
//        sim.integrator.getEventManager().addListener(configStorage);
//        configStorage.addListener(meterBAC);
//        configStorage.addListener(meterRAC);
//        configStorage.addListener(meterMSDAC);

        AccumulatorAverageCovariance accumulatorCentVir = new AccumulatorAverageCovariance(blockSize);
        DataPumpListener accumulatorPumpCentVir = new DataPumpListener(meterCentVir, accumulatorCentVir, interval);
        sim.integrator.getEventManager().addListener(accumulatorPumpCentVir);

        //Short run for Cv calculations, to avoid <Big^2>-<Big>^2 roundoff.
        double EnShift = 0, errEnShift = 0;
        long numStepsShort = steps/10;
        System.out.println(" Short sim for Covariance: " + numStepsShort + " numSteps");
        sim.getController().runActivityBlocking(new ActivityIntegrate(sim.integrator, numStepsShort));
        System.out.println(" Done with "+ numStepsShort +" steps of short run");
        IData dataAvgShort = accumulatorCentVir.getData(accumulatorCentVir.AVERAGE);
        IData dataErrShort = accumulatorCentVir.getData(accumulatorCentVir.ERROR);
        EnShift = dataAvgShort.getValue(0);
        errEnShift = dataErrShort.getValue(0);
        System.out.println(" EnShift:  " + EnShift/numAtoms + "  err: " + errEnShift/numAtoms);
        accumulatorCentVir.reset();
        meterCentVir.setEnShift(EnShift);


        AccumulatorAverageCovariance accumulatorKE = new AccumulatorAverageCovariance(blockSize);
        DataPumpListener accumulatorPumpKE = new DataPumpListener(meterKE, accumulatorKE, interval);
        sim.integrator.getEventManager().addListener(accumulatorPumpKE);

        // Primitive
        AccumulatorAverageCovariance accumulatorPrim = new AccumulatorAverageCovariance(blockSize);
        if (meterPrim != null) {
            meterPrim.setEnShift(EnShift);
            DataPumpListener accumulatorPumpPrim = new DataPumpListener(meterPrim, accumulatorPrim, interval);
            sim.integrator.getEventManager().addListener(accumulatorPumpPrim);
        }

        // Virial
        AccumulatorAverageCovariance accumulatorVir = new AccumulatorAverageCovariance(blockSize);
        if (meterVir != null) {
            meterVir.setEnShift(EnShift);
            DataPumpListener accumulatorPumpVir = new DataPumpListener(meterVir, accumulatorVir, interval);
            sim.integrator.getEventManager().addListener(accumulatorPumpVir);
        }

        AccumulatorAverageCovariance accumulatorCentVirFD = new AccumulatorAverageCovariance(blockSize);
        if (meterCentVirFD != null) {
            meterCentVirFD.setEnShift(EnShift);
            DataPumpListener accumulatorPumpCentVirFD = new DataPumpListener(meterCentVirFD, accumulatorCentVirFD, interval);
            sim.integrator.getEventManager().addListener(accumulatorPumpCentVirFD);
        }

        // HMAc
        AccumulatorAverageCovariance accumulatorHMAc = new AccumulatorAverageCovariance(blockSize);
        if (meterHMAc != null) {
            meterHMAc.setEnShift(EnShift);
            DataPumpListener accumulatorPumpHMAc = new DataPumpListener(meterHMAc, accumulatorHMAc, interval);
            sim.integrator.getEventManager().addListener(accumulatorPumpHMAc);
        }

        // HMAc FD
        AccumulatorAverageCovariance accumulatorHMAcFD = new AccumulatorAverageCovariance(blockSize);
        if (meterHMAcFD != null) {
            meterHMAcFD.setEnShift(EnShift);
            DataPumpListener accumulatorPumpHMAcFD = new DataPumpListener(meterHMAcFD, accumulatorHMAcFD, interval);
            sim.integrator.getEventManager().addListener(accumulatorPumpHMAcFD);
        }

        //Simple NM
        AccumulatorAverageCovariance accumulatorNMSimple = new AccumulatorAverageCovariance(blockSize);
        if (meterNMSimple != null) {
            meterNMSimple.setEnShift(EnShift);
            DataPumpListener accumulatorPumpHMAsimple = new DataPumpListener(meterNMSimple, accumulatorNMSimple, interval);
            sim.integrator.getEventManager().addListener(accumulatorPumpHMAsimple);
        }

        //HO NM
        AccumulatorAverageCovariance accumulatorNMEC = new AccumulatorAverageCovariance(blockSize);
        if (meterNMEC != null) {
            meterNMEC.setEnShift(EnShift);
            DataPumpListener accumulatorPumpNMEC = new DataPumpListener(meterNMEC, accumulatorNMEC, interval);
            sim.integrator.getEventManager().addListener(accumulatorPumpNMEC);
        }

        //Simple stage
        AccumulatorAverageCovariance accumulatorStageSimple = new AccumulatorAverageCovariance(blockSize);
        if (meterStageSimple != null) {
            meterStageSimple.setEnShift(EnShift);
            DataPumpListener pumpStagesimple = new DataPumpListener(meterStageSimple, accumulatorStageSimple, interval);
            sim.integrator.getEventManager().addListener(pumpStagesimple);
        }

        //EC stage
        AccumulatorAverageCovariance accumulatorStageEC = new AccumulatorAverageCovariance(blockSize);
        if (meterStageEC != null) {
            meterStageEC.setEnShift(EnShift);
            DataPumpListener pumpStageEC = new DataPumpListener(meterStageEC, accumulatorStageEC, interval);
            sim.integrator.getEventManager().addListener(pumpStageEC);
        }

        //Run ...
        sim.integrator.resetStepCount();
        sim.getController().runActivityBlocking(new ActivityIntegrate(sim.integrator, steps));

        //T
        DataGroup dataKE = (DataGroup) accumulatorKE.getData();
        double avgKE = dataKE.getValue(accumulatorKE.AVERAGE.index);
        double errKE = dataKE.getValue(accumulatorKE.ERROR.index);
        double corKE = dataKE.getValue(accumulatorKE.BLOCK_CORRELATION.index);
        System.out.println(" T_measured: " + avgKE / (1.5 * numAtoms * nBeads) + " +/- " + errKE / (1.5 * numAtoms * nBeads) + " cor: " + corKE);
        System.out.println();
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
//        double avgPnCentVir = dataAvgCentVir.getValue(2);
//        double errPnCentVir = dataErrCentVir.getValue(2);
//        double corPnCentVir = dataCorCentVir.getValue(2);
        double CvnCentVir = kB_beta2*(dataAvgCentVir.getValue(1) - avgEnCentVir*avgEnCentVir);
        varX0 = errEnCentVir*errEnCentVir;
        varX1 = dataErrCentVir.getValue(1)*dataErrCentVir.getValue(1);
        corX0X1 = dataCovCentVir.getValue(1)/Math.sqrt(dataCovCentVir.getValue(0))/Math.sqrt(dataCovCentVir.getValue(3));
        double errCvnCentVir = kB_beta2*Math.sqrt(varX1 + 4.0*avgEnCentVir*avgEnCentVir*varX0 - 4*avgEnCentVir*dataErrCentVir.getValue(0)*dataErrCentVir.getValue(1)*corX0X1);

        // Cent-Vir: FD
        DataGroup dataCentVirFD = (DataGroup) accumulatorCentVirFD.getData();
        IData dataAvgCentVirFD = dataCentVirFD.getData(accumulatorCentVirFD.AVERAGE.index);
        IData dataErrCentVirFD = dataCentVirFD.getData(accumulatorCentVirFD.ERROR.index);
        IData dataCorCentVirFD = dataCentVirFD.getData(accumulatorCentVirFD.BLOCK_CORRELATION.index);
        IData dataCovCentVirFD = dataCentVirFD.getData(accumulatorCentVirFD.COVARIANCE.index);
        double avgEnCentVirFD = dataAvgCentVirFD.getValue(0);
        double errEnCentVirFD = dataErrCentVirFD.getValue(0);
        double corEnCentVirFD = dataCorCentVirFD.getValue(0);
//        double avgPnCentVirFD = dataAvgCentVirFD.getValue(2);
//        double errPnCentVirFD = dataErrCentVirFD.getValue(2);
//        double corPnCentVirFD = dataCorCentVirFD.getValue(2);
        double CvnCentVirFD = kB_beta2*(dataAvgCentVirFD.getValue(1) - avgEnCentVirFD*avgEnCentVirFD);
        varX0 = errEnCentVirFD*errEnCentVirFD;
        varX1 = dataErrCentVirFD.getValue(1)*dataErrCentVirFD.getValue(1);
        corX0X1 = dataCovCentVirFD.getValue(1)/Math.sqrt(dataCovCentVirFD.getValue(0))/Math.sqrt(dataCovCentVirFD.getValue(3));
        double errCvnCentVirFD = kB_beta2*Math.sqrt(varX1 + 4.0*avgEnCentVirFD*avgEnCentVirFD*varX0 - 4*avgEnCentVirFD*dataErrCentVirFD.getValue(0)*dataErrCentVirFD.getValue(1)*corX0X1);

        double avgEnPrim=0, avgEnVir=0, avgEnHMAc=0, avgEnNMSimple=0, avgEnNMEC=0, avgEnStageSimple=0,avgEnStageEC=0;
        double errEnPrim=0, errEnVir=0, errEnHMAc=0, errEnNMSimple=0, errEnNMEC=0, errEnStageSimple=0, errEnStageEC=0;
        double corEnPrim=0, corEnVir=0, corEnHMAc=0, corEnNMSimple=0, corEnNMEC=0, corEnStageSimple=0, corEnStageEC=0;
        double CvnPrim=0, CvnVir=0, CvnHMAc=0, CvnNMSimple=0, CvnNMEC=0, CvnStageSimple=0, CvnStageEC=0;
        double errCvnPrim=0, errCvnVir=0, errCvnHMAc=0, errCvnNMsimple=0, errCvnNMEC=0, errCvnStageSimple=0, errCvnStageEC=0;
        double avgPnPrim = 0, errPnPrim = 0, corPnPrim = 0;

        DataGroup dataPrim = (DataGroup) accumulatorPrim.getData();
        IData dataAvgPrim = dataPrim.getData(accumulatorPrim.AVERAGE.index);
        IData dataErrPrim = dataPrim.getData(accumulatorPrim.ERROR.index);
        IData dataCorPrim = dataPrim.getData(accumulatorPrim.BLOCK_CORRELATION.index);
        IData dataCovPrim = dataPrim.getData(accumulatorPrim.COVARIANCE.index);
        avgEnPrim = dataAvgPrim.getValue(0);
        errEnPrim = dataErrPrim.getValue(0);
        corEnPrim = dataCorPrim.getValue(0);
//        avgPnPrim = dataAvgPrim.getValue(2);
//        errPnPrim = dataErrPrim.getValue(2);
//        corPnPrim = dataCorPrim.getValue(2);
        CvnPrim = kB_beta2*(dataAvgPrim.getValue(1) - avgEnPrim*avgEnPrim);
        varX0 = errEnPrim*errEnPrim;
        varX1 = dataErrPrim.getValue(1)*dataErrPrim.getValue(1);
        corX0X1 = dataCovPrim.getValue(1)/Math.sqrt(dataCovPrim.getValue(0))/Math.sqrt(dataCovPrim.getValue(3));
        errCvnPrim = kB_beta2*Math.sqrt(varX1 + 4.0*avgEnPrim*avgEnPrim*varX0 - 4*avgEnPrim*dataErrPrim.getValue(0)*dataErrPrim.getValue(1)*corX0X1);


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


        // HMAc: FD
        DataGroup dataHMAcFD = (DataGroup) accumulatorHMAcFD.getData();
        IData dataAvgHMAcFD = dataHMAcFD.getData(accumulatorHMAcFD.AVERAGE.index);
        IData dataErrHMAcFD = dataHMAcFD.getData(accumulatorHMAcFD.ERROR.index);
        IData dataCorHMAcFD = dataHMAcFD.getData(accumulatorHMAcFD.BLOCK_CORRELATION.index);
        IData dataCovHMAcFD = dataHMAcFD.getData(accumulatorHMAcFD.COVARIANCE.index);
        double avgEnHMAcFD = dataAvgHMAcFD.getValue(0) ;
        double errEnHMAcFD = dataErrHMAcFD.getValue(0);
        double corEnHMAcFD = dataCorHMAcFD.getValue(0);
        double CvnHMAcFD = kB_beta2*(dataAvgHMAcFD.getValue(1) - avgEnHMAcFD*avgEnHMAcFD);
        varX0 = errEnHMAcFD*errEnHMAcFD;
        varX1 = dataErrHMAcFD.getValue(1)*dataErrHMAcFD.getValue(1);
        corX0X1 = dataCovHMAcFD.getValue(1)/Math.sqrt(dataCovHMAcFD.getValue(0))/Math.sqrt(dataCovHMAcFD.getValue(3));
        double errCvnHMAcFD = kB_beta2*Math.sqrt(varX1 + 4.0*avgEnHMAcFD*avgEnHMAcFD*varX0 - 4*avgEnHMAcFD*dataErrHMAcFD.getValue(0)*dataErrHMAcFD.getValue(1)*corX0X1);
        if (errEnHMAcFD < 1e-10){
            errCvnHMAcFD = kB_beta2*Math.sqrt(varX1 + 4.0*avgEnHMAcFD*avgEnHMAcFD*varX0 );
        }


        if (!onlyCentroid) {
            // Prim
//            DataGroup dataPrim = (DataGroup) accumulatorPrim.getData();
//            IData dataAvgPrim = dataPrim.getData(accumulatorPrim.AVERAGE.index);
//            IData dataErrPrim = dataPrim.getData(accumulatorPrim.ERROR.index);
//            IData dataCorPrim = dataPrim.getData(accumulatorPrim.BLOCK_CORRELATION.index);
//            IData dataCovPrim = dataPrim.getData(accumulatorPrim.COVARIANCE.index);
//            avgEnPrim = dataAvgPrim.getValue(0);
//            errEnPrim = dataErrPrim.getValue(0);
//            corEnPrim = dataCorPrim.getValue(0);
//            CvnPrim = kB_beta2*(dataAvgPrim.getValue(1) - avgEnPrim*avgEnPrim);
//            varX0 = errEnPrim*errEnPrim;
//            varX1 = dataErrPrim.getValue(1)*dataErrPrim.getValue(1);
//            corX0X1 = dataCovPrim.getValue(1)/Math.sqrt(dataCovPrim.getValue(0))/Math.sqrt(dataCovPrim.getValue(3));
//            errCvnPrim = kB_beta2*Math.sqrt(varX1 + 4.0*avgEnPrim*avgEnPrim*varX0 - 4*avgEnPrim*dataErrPrim.getValue(0)*dataErrPrim.getValue(1)*corX0X1);

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
            System.out.println(" En_prim:         " + (avgEnPrim+EnShift)/numAtoms + "   err: " + errEnPrim/numAtoms + " cor: " + corEnPrim);
            System.out.println(" En_cvir:         " + (avgEnCentVir+EnShift)/numAtoms + "   err: " + errEnCentVir/numAtoms + " cor: " + corEnCentVir);
            System.out.println(" En_cvir_fd:      " + (avgEnCentVirFD+EnShift)/numAtoms + "   err: " + errEnCentVirFD/numAtoms + " cor: " + corEnCentVirFD);
            System.out.println(" En_hmac:         " + (avgEnHMAc+EnShift)/numAtoms + "   err: " + errEnHMAc/numAtoms + " cor: " + corEnHMAc);
            System.out.println(" En_hmac_fd:      " + (avgEnHMAcFD+EnShift)/numAtoms + "   err: " + errEnHMAcFD/numAtoms + " cor: " + corEnHMAcFD);
            System.out.println();
//            System.out.println(" Pn_prim:         " + avgPnPrim + "   err: " + errPnPrim + " cor: " + corPnPrim);
//            System.out.println(" Pn_cvir:         " + avgPnCentVir + "   err: " + errPnCentVir + " cor: " + corPnCentVir);
//            System.out.println();
            System.out.println(" Cvn_prim:         " + CvnPrim/numAtoms +       "   err: " + errCvnPrim/numAtoms);
            System.out.println(" Cvn_cvir:         " + CvnCentVir/numAtoms +    "   err: " + errCvnCentVir/numAtoms);
            System.out.println(" Cvn_cvir_fd:      " + CvnCentVirFD/numAtoms +  "   err: " + errCvnCentVirFD/numAtoms);
            System.out.println(" Cvn_hmac:         " + CvnHMAc/numAtoms +       "   err: " + errCvnHMAc/numAtoms);
            System.out.println(" Cvn_hmac_fd:      " + CvnHMAcFD/numAtoms +     "   err: " + errCvnHMAcFD/numAtoms);

        } else {
            System.out.println("\n En_prim:          " + avgEnPrim + "   err: " + errEnPrim + " cor: " + corEnPrim);
            System.out.println(" En_vir:          " + (avgEnVir+EnShift)/numAtoms + "   err: " + errEnVir/numAtoms + " cor: " + corEnVir);
            System.out.println(" En_nm_simple:    " + (avgEnNMSimple+EnShift)/numAtoms + "   err: " + errEnNMSimple/numAtoms + " cor: " + corEnNMSimple);
            System.out.println(" En_nm_ec:        " + (avgEnNMEC+EnShift)/numAtoms + "   err: " + errEnNMEC/numAtoms + " cor: " + corEnNMEC);
            System.out.println(" En_stage_simple: " + (avgEnStageSimple+EnShift)/numAtoms + "   err: " + errEnStageSimple/numAtoms + " cor: " + corEnStageSimple);
            System.out.println(" En_stage_ec:     " + (avgEnStageEC+EnShift)/numAtoms + "   err: " + errEnStageEC/numAtoms + " cor: " + corEnStageEC);
            System.out.println();
            System.out.println(" Cvn_prim:         " + CvnPrim/numAtoms +          "   err: " + errCvnPrim/numAtoms);
            System.out.println(" Cvn_nm_simple:    " + CvnNMSimple/numAtoms +    "   err: " + errCvnNMsimple/numAtoms);
            System.out.println(" Cvn_nm_ec:        " + CvnNMEC/numAtoms +          "   err: " + errCvnNMEC/numAtoms);
            System.out.println(" Cvn_stage_simple: " + CvnStageSimple/numAtoms + "   err: " + errCvnStageSimple/numAtoms);
            System.out.println(" Cvn_stage_ec:     " + CvnStageEC/numAtoms +       "   err: " + errCvnStageEC/numAtoms);
        }


        long t2 = System.currentTimeMillis();
//        try {
//            writeDataToFile(meterBAC, meterBAC.getErrorData(), "bac.dat");
//            writeDataToFile(meterRAC, meterRAC.getErrorData(), "rac.dat");
//            writeDataToFile(meterMSDAC, meterMSDAC.getErrorData(), "msdac.dat");
//        } catch (IOException e) {
//            throw new RuntimeException(e);
//        }

        System.out.println("\n time: (min) " + (t2 - t1) * 0.001 / 60.0);
    }

    public static void writeDataToFile(IDataSource meter, IData errData, String filename) throws IOException {
        IData data;
        IData xData;
        if (meter instanceof AccumulatorAverage) {
            AccumulatorAverage acc = (AccumulatorAverage) meter;
            data = acc.getData(acc.AVERAGE);
            xData = ((DataFunction.DataInfoFunction) ((DataGroup.DataInfoGroup) acc.getDataInfo()).getSubDataInfo(acc.AVERAGE.index)).getXDataSource().getIndependentData(0);
        } else {
            data = meter.getData();
            xData = ((DataFunction.DataInfoFunction) meter.getDataInfo()).getXDataSource().getIndependentData(0);
        }
        boolean allNaN = true;
        for (int i = 0; i < xData.getLength(); i++) {
            if (!Double.isNaN(data.getValue(i))) allNaN = false;
        }
        if (allNaN) return;
        FileWriter fw = new FileWriter(filename);
        for (int i = 0; i < xData.getLength(); i++) {
            double y = data.getValue(i);
            if (Double.isNaN(y)) continue;
            if (errData == null) {
                fw.write(xData.getValue(i) + " " + y + "\n");
            } else {
                fw.write(xData.getValue(i) + " " + y + " " + errData.getValue(i) + "\n");
            }
        }
        fw.close();
    }

    public enum MoveChoice {Real, NM, NMEC, Stage, StageEC};

    public static class SimParams extends ParameterBase {
        public int D = 3;
        public double k2 = 219.949835;
        public double gammaLangevin = Math.sqrt(k2);
        public long steps = 100000;
        public double density = 1.0;
        public double temperature = 1.0;
        public int numAtoms = 108;
        public double mass = 1.0;
        public double hbar = 0.1;
        public double rc = 2.5;
        public boolean isGraphic = false;
        public int nShifts = 0;
        public double timeStep = -1;
        public int nBeads = -1;
        public MoveChoice coordType = MoveChoice.StageEC;
        public double facTimestep = 0.01;
        public boolean onlyCentroid = true;
    }
}