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
import etomica.data.meter.MeterEnergyFromIntegrator;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.data.meter.MeterPotentialEnergyFromIntegrator;
import etomica.data.types.DataFunction;
import etomica.data.types.DataGroup;
import etomica.graphics.*;
import etomica.integrator.IntegratorLangevin;
import etomica.integrator.IntegratorListener;
import etomica.integrator.IntegratorMD;
import etomica.lattice.LatticeCubicFcc;
import etomica.nbr.list.PotentialMasterList;
import etomica.potential.*;
import etomica.potential.compute.PotentialCompute;
import etomica.potential.compute.PotentialComputeAggregate;
import etomica.potential.compute.PotentialComputePair;
import etomica.simulation.Simulation;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.space3d.Space3D;
import etomica.species.SpeciesBuilder;
import etomica.species.SpeciesGeneral;
import etomica.units.dimensions.Length;
import etomica.util.ParameterBase;
import etomica.util.ParseArgs;
import java.awt.*;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import etomica.util.random.RandomMersenneTwister;
import org.apache.batik.css.engine.value.svg.StrokeWidthManager;

/**
 * PIMD using Langevin Thermostat using BAOAB algorithm.
 *
 * @author Andrew Schultz
 */

public class LJPIMD extends Simulation {

    public final PotentialComputePair potentialMaster;
    public final IntegratorMD integrator;
    public final PotentialMasterBonding pmBonding;
    public final Box box;
    public final PotentialComputeAggregate pmAgg;
    public final MCMoveHOReal2 moveStageSimple, moveStageEC;
    public double beta;
    public int dim;

    public LJPIMD(Space space, MoveChoice coordType, double mass, double timeStep, double gammaLangevin, int nBeads, int numAtoms, double temperature, double density, double rc, double omega2, double hbar) {
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
        NeighborListManagerPI neighborManager = new NeighborListManagerPI(getSpeciesManager(), box, 2, rc, BondingInfo.noBonding());
        neighborManager.setAutoUpdateNeighbors(false);
        potentialMaster = new PotentialComputePair(getSpeciesManager(), box, neighborManager);
        pmBonding = new PotentialMasterBonding(getSpeciesManager(), box);
        beta = 1.0/temperature;

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
        P2SoftSphericalTruncated p2 = new P2SoftSphericalTruncated(p2lj, rc);
        AtomType atomType = species.getLeafType();
        potentialMaster.setPairPotential(atomType, atomType, p2);
        potentialMaster.doAllTruncationCorrection = false;

        PotentialComputeAggregate.localStorageDefault = true;
        pmAgg = new PotentialComputeAggregate(pmBonding, potentialMaster);

        if (coordType == LJPIMD.MoveChoice.Real) {
            integrator = new IntegratorLangevin(pmAgg, random, timeStep, temperature, box, gammaLangevin);
        } else if (coordType == LJPIMD.MoveChoice.NM) {
            MCMoveHO move = new MCMoveHO(space, pmAgg, random, temperature, 0, box, hbar);
            integrator = new IntegratorLangevinPINM(pmAgg, random, timeStep, temperature, box, gammaLangevin, move, hbar, omega2);
        } else if (coordType == LJPIMD.MoveChoice.NMEC) {
            MCMoveHO move = new MCMoveHO(space, pmAgg, random, temperature, omega2, box, hbar);
            integrator = new IntegratorLangevinPINM(pmAgg, random, timeStep, temperature, box, gammaLangevin, move, hbar, omega2);
        } else if (coordType == LJPIMD.MoveChoice.Stage) {
            MCMoveHOReal2 move = new MCMoveHOReal2(space, pmAgg, random, temperature, 0, box, hbar);
            integrator = new IntegratorLangevinPI(pmAgg, random, timeStep, temperature, box, gammaLangevin, move, hbar, omega2);
        } else { //StageEC -- default
            MCMoveHOReal2 move = new MCMoveHOReal2(space, pmAgg, random, temperature, omega2, box, hbar);
            integrator = new IntegratorLangevinPI(pmAgg, random, timeStep, temperature, box, gammaLangevin, move, hbar, omega2);
        }

        // find neighbors now.  Don't hook up NeighborListManager (neighbors won't change)
        potentialMaster.init();
        // and then never update neighbor lists
        integrator.getEventManager().removeListener(neighborManager);
        p2.setTruncationRadius(2*rc);

        moveStageSimple = new MCMoveHOReal2(space, pmAgg, random, temperature, 0, box, hbar);
        moveStageEC = new MCMoveHOReal2(space, pmAgg, random, temperature, omega2, box, hbar);
        integrator.setThermostatNoDrift(false);
        integrator.setIsothermal(true);
    }

    public static void main(String[] args) {
        long t1 = System.currentTimeMillis();
        LJPIMD.SimParams params = new LJPIMD.SimParams();
        if (args.length > 0) {
            ParseArgs.doParseArgs(params, args);
        } else {
            params.steps = 10000;
            params.hbar = 0.1;
            params.temperature = 0.5;
            params.numAtoms = 32;
            params.rc = 3;
            params.isGraphic = false;
            params.onlyCentroid = false;
//            params.coordType = MoveChoice.Real;
//            params.coordType = MoveChoice.NM;
//            params.coordType = MoveChoice.Stage;
//            params.coordType = MoveChoice.NMEC;
            params.coordType = LJPIMD.MoveChoice.StageEC;
            params.timeStep = 0.01;

            params.nBeads = 1;


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
        LJPIMD.MoveChoice coordType = params.coordType;
        long steps = params.steps;
        long stepsEq = steps/10;

        double omega = Math.sqrt(omega2);
        double x = 1/temperature*hbar*omega;
        if (nBeads == -1){
            nBeads = (int) (20*x);
        }
        double omegaN = nBeads*temperature/hbar;

        LJPIMD sim = new LJPIMD(space, coordType, mass, timeStep, gammaLangevin, nBeads, numAtoms, temperature, density, rc, omega2, hbar);
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
        MeterPICentVir meterCentVir = null;
        MeterPIHMAc meterHMAc = null;
        MeterPIHMA meterNMEC = null;
        MeterPIHMAReal2 meterStageEC = null;

        meterPrim = new MeterPIPrim(sim.pmBonding, sim.potentialMaster, temperature, nBeads, sim.box);
        meterCentVir = new MeterPICentVir(sim.potentialMaster, temperature, nBeads, sim.box);
        meterHMAc = new MeterPIHMAc(sim.potentialMaster, temperature, nBeads, sim.box);
        meterNMEC = new MeterPIHMA(sim.pmBonding, sim.potentialMaster, temperature, nBeads, omega2, sim.box, hbar);
        meterStageEC = new MeterPIHMAReal2(sim.pmBonding, sim.potentialMaster, temperature, nBeads, sim.moveStageEC, nShifts);


        int interval = 5;
        int blocks = 100;
        long blockSize = steps / (interval * blocks);
        if (blockSize == 0) blockSize = 1;
        System.out.println(" numBlocks: " + blocks + " blocksize: " + blockSize + " interval: " + interval);

        if (isGraphic) {
            sim.getController().addActivity(new ActivityIntegrate(sim.integrator));
            SimulationGraphic simGraphic = new SimulationGraphic(sim, SimulationGraphic.GRAPHIC_ONLY);
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
            ((DiameterHashByType) ((DisplayBox) simGraphic.displayList().getFirst()).getDiameterHash()).setDiameter(sim.species().getAtomType(0), 0.1);
            simGraphic.makeAndDisplayFrame("PIMD-"+coordType);

            return;
        }

        System.out.flush();

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

        // HMAc
        AccumulatorAverageCovariance accumulatorHMAc = new AccumulatorAverageCovariance(blockSize);
        if (meterHMAc != null) {
            meterHMAc.setEnShift(EnShift);
            DataPumpListener accumulatorPumpHMAc = new DataPumpListener(meterHMAc, accumulatorHMAc, interval);
            sim.integrator.getEventManager().addListener(accumulatorPumpHMAc);
        }

        //NMEC
        AccumulatorAverageCovariance accumulatorNMEC = new AccumulatorAverageCovariance(blockSize);
        if (meterNMEC != null) {
            meterNMEC.setEnShift(EnShift);
            DataPumpListener accumulatorPumpNMEC = new DataPumpListener(meterNMEC, accumulatorNMEC, interval);
            sim.integrator.getEventManager().addListener(accumulatorPumpNMEC);
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
        double kB_beta2 = sim.beta*sim.beta;
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

        double avgEnPrim=0, avgEnVir=0, avgEnHMAc=0, avgEnNMEC=0;
        double errEnPrim=0, errEnVir=0, errEnHMAc=0, errEnNMEC=0;
        double corEnPrim=0, corEnVir=0, corEnHMAc=0, corEnNMEC=0;
        double CvnPrim=0, CvnVir=0, CvnHMAc=0, CvnNMEC=0;
        double errCvnPrim=0, errCvnVir=0, errCvnHMAc=0, errCvnNMEC=0;
        double avgPnPrim = 0, errPnPrim = 0, corPnPrim = 0;

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

        // NMEC
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

        // HMA EC stage
        DataGroup dataStageEC = (DataGroup) accumulatorStageEC.getData();
        IData dataAvgStageEC = dataStageEC.getData(accumulatorStageEC.AVERAGE.index);
        IData dataErrStageEC = dataStageEC.getData(accumulatorStageEC.ERROR.index);
        IData dataCorStageEC = dataStageEC.getData(accumulatorStageEC.BLOCK_CORRELATION.index);
        IData dataCovStageEC = dataStageEC.getData(accumulatorStageEC.COVARIANCE.index);
        double avgEnStageEC = dataAvgStageEC.getValue(0);
        double errEnStageEC = dataErrStageEC.getValue(0);
        double corEnStageEC = dataCorStageEC.getValue(0);
        double CvnStageEC  = kB_beta2*(dataAvgStageEC.getValue(1) - avgEnStageEC*avgEnStageEC);
        varX0 = errEnStageEC*errEnStageEC;
        varX1 = dataErrStageEC.getValue(1)*dataErrStageEC.getValue(1);
        corX0X1 = dataCovStageEC.getValue(1)/Math.sqrt(dataCovStageEC.getValue(0))/Math.sqrt(dataCovStageEC.getValue(3));
        double errCvnStageEC = kB_beta2*Math.sqrt(varX1 + 4.0*avgEnStageEC*avgEnStageEC*varX0 - 4*avgEnStageEC*errEnStageEC*dataErrStageEC.getValue(1)*corX0X1);
        if (errEnStageEC < 1e-10){
            errCvnStageEC = kB_beta2*Math.sqrt(varX1 + 4.0*avgEnStageEC*avgEnStageEC*varX0);
        }

        System.out.println(" En_prim:        " + (avgEnPrim+EnShift)/numAtoms + "   err: " + errEnPrim/numAtoms + " cor: " + corEnPrim);
        System.out.println(" En_cvir:        " + (avgEnCentVir+EnShift)/numAtoms + "   err: " + errEnCentVir/numAtoms + " cor: " + corEnCentVir);
        System.out.println(" En_hmac:        " + (avgEnHMAc+EnShift)/numAtoms + "   err: " + errEnHMAc/numAtoms + " cor: " + corEnHMAc);
        System.out.println(" En_nm_ec:       " + (avgEnNMEC+EnShift)/numAtoms + "   err: " + errEnNMEC/numAtoms + " cor: " + corEnNMEC);
        System.out.println(" En_stage_ec:    " + (avgEnStageEC+EnShift)/numAtoms + "   err: " + errEnStageEC/numAtoms + " cor: " + corEnStageEC);

        System.out.println();
        System.out.println(" Cvn_prim:         " + CvnPrim/numAtoms +       "   err: " + errCvnPrim/numAtoms);
        System.out.println(" Cvn_cvir:         " + CvnCentVir/numAtoms +    "   err: " + errCvnCentVir/numAtoms);
        System.out.println(" Cvn_hmac:         " + CvnHMAc/numAtoms +       "   err: " + errCvnHMAc/numAtoms);
        System.out.println(" Cvn_nm_ec:        " + CvnNMEC/numAtoms +       "   err: " + errCvnNMEC/numAtoms);
        System.out.println(" Cvn_stage_ec:        " + CvnStageEC/numAtoms + "   err: " + errCvnStageEC/numAtoms);

        long t2 = System.currentTimeMillis();
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
        public double k2 = 218.220183;
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
        public LJPIMD.MoveChoice coordType = LJPIMD.MoveChoice.StageEC;
        public double facTimestep = 0.01;
        public boolean onlyCentroid = true;
    }
}
