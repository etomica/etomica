/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.normalmode;

import etomica.action.activity.ActivityIntegrate;
import etomica.atom.AtomType;
import etomica.box.Box;
import etomica.data.*;
import etomica.data.histogram.HistogramSimple;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.data.meter.MeterPotentialEnergyFromIntegrator;
import etomica.data.types.DataDoubleArray;
import etomica.data.types.DataGroup;
import etomica.integrator.IntegratorBox;
import etomica.integrator.IntegratorMC;
import etomica.integrator.mcmove.MCMoveStepTracker;
import etomica.lattice.crystal.*;
import etomica.math.DoubleRange;
import etomica.overlap.IntegratorOverlap;
import etomica.potential.P2SoftSphere;
import etomica.potential.P2SoftSphericalTruncatedShifted;
import etomica.potential.Potential2SoftSpherical;
import etomica.potential.PotentialMasterMonatomic;
import etomica.simulation.Simulation;
import etomica.space.Boundary;
import etomica.space.BoundaryRectangularPeriodic;
import etomica.space.Space;
import etomica.species.SpeciesSpheresMono;
import etomica.util.ParameterBase;
import etomica.util.ReadParameters;
import etomica.virial.overlap.AccumulatorVirialOverlapSingleAverage;
import etomica.virial.overlap.DataSourceVirialOverlap;

import java.io.*;

/**
 * Simulation to run sampling with the hard sphere potential, but measuring
 * the harmonic potential based on normal mode data from a previous simulation.
 * 
 * To study the probability distribution of P_aw/ P_wa and P_bw/ P_wb
 * a, b and w being the harmonic, target and the overlap region respectively
 * 
 * From there, we compute the relative entropy term.
 * 
 * @author  Tai Boon Tan
 */
public class SimOverlapSoftSphereReweighting extends Simulation {

    private static final long serialVersionUID = 1L;
    public final double latticeEnergy;
    public IntegratorOverlap integratorOverlap;
    public DataSourceVirialOverlap dsvo;
    public IntegratorBox[] integrators;
    public ActivityIntegrate activityIntegrate;
    public Box boxTarget, boxHarmonic;
    public Boundary boundaryTarget, boundaryHarmonic;
    public int[] nCells;
    public Basis basis;
    public NormalModes normalModes;
    public Primitive primitive;
    public double refPref;
    public AccumulatorVirialOverlapSingleAverage[] accumulators;
    public MeterHarmonicEnergy meterHarmonicEnergy;
    public DataPumpListener[] accumulatorPumps;
    public IDataSource[] meters;
    public String fname;
    protected MCMoveHarmonic move;
    protected PotentialMasterMonatomic potentialMasterTarget;
    public SimOverlapSoftSphereReweighting(Space _space, int numAtoms, double density, double temperature, String filename, double harmonicFudge, int exponent) {
        super(_space);
        this.fname = filename;
        potentialMasterTarget = new PotentialMasterMonatomic(this);
        integrators = new IntegratorBox[2];
        accumulatorPumps = new DataPumpListener[2];
        meters = new IDataSource[2];
        accumulators = new AccumulatorVirialOverlapSingleAverage[2];

        SpeciesSpheresMono species = new SpeciesSpheresMono(this, space);
        addSpecies(species);

        // TARGET

        boxTarget = this.makeBox();
        boxTarget.setNMolecules(species, numAtoms);

        IntegratorMC integratorTarget = new IntegratorMC(potentialMasterTarget, getRandom(), temperature, boxTarget);
        MCMoveAtomCoupled atomMove = new MCMoveAtomCoupled(potentialMasterTarget, new MeterPotentialEnergy(potentialMasterTarget), getRandom(), space);
        atomMove.setStepSize(0.1);
        atomMove.setStepSizeMax(0.5);
        integratorTarget.getMoveManager().addMCMove(atomMove);
        ((MCMoveStepTracker) atomMove.getTracker()).setNoisyAdjustment(true);

        integrators[1] = integratorTarget;

        if (space.D() == 1) {
            primitive = new PrimitiveCubic(space, 1.0 / density);
            boundaryTarget = new BoundaryRectangularPeriodic(space, numAtoms / density);
            nCells = new int[]{numAtoms};
            basis = new BasisMonatomic(space);
        } else {
            double L = Math.pow(4.0 / density, 1.0 / 3.0);
            primitive = new PrimitiveCubic(space, L);
            int n = (int) Math.round(Math.pow(numAtoms / 4, 1.0 / 3.0));
            nCells = new int[]{n, n, n};
            boundaryTarget = new BoundaryRectangularPeriodic(space, n * L);
            basis = new BasisCubicFcc();
        }
        boxTarget.setBoundary(boundaryTarget);

        CoordinateDefinitionLeaf coordinateDefinitionTarget = new CoordinateDefinitionLeaf(boxTarget, primitive, basis, space);
        coordinateDefinitionTarget.initializeCoordinates(nCells);

        Potential2SoftSpherical potential = new P2SoftSphere(space, 1.0, 1.0, exponent);
        double truncationRadius = boundaryTarget.getBoxSize().getX(0) * 0.495;
        P2SoftSphericalTruncatedShifted pTruncated = new P2SoftSphericalTruncatedShifted(space, potential, truncationRadius);
        AtomType sphereType = species.getLeafType();
        potentialMasterTarget.addPotential(pTruncated, new AtomType[]{sphereType, sphereType});
        atomMove.setPotential(pTruncated);

        /*
         *  1-body Potential to Constraint the atom from moving too far
         *  	away from its lattice-site
         *
         */

        P1Constraint p1Constraint = new P1Constraint(space, primitive.getSize()[0], boxTarget, coordinateDefinitionTarget);
        potentialMasterTarget.addPotential(p1Constraint, new AtomType[]{sphereType});

        potentialMasterTarget.lrcMaster().setEnabled(false);
        MeterPotentialEnergy meterPE = new MeterPotentialEnergy(potentialMasterTarget, boxTarget);
        latticeEnergy = meterPE.getDataAsScalar();


        // HARMONIC
        boundaryHarmonic = new BoundaryRectangularPeriodic(space);
        boxHarmonic = this.makeBox(boundaryHarmonic);
        boxHarmonic.setNMolecules(species, numAtoms);

        IntegratorMC integratorHarmonic = new IntegratorMC(potentialMasterTarget, random, 1.0, boxHarmonic);

        move = new MCMoveHarmonic(getRandom());
        integratorHarmonic.getMoveManager().addMCMove(move);
        integrators[0] = integratorHarmonic;

        if (space.D() == 1) {
            boundaryHarmonic = new BoundaryRectangularPeriodic(space, numAtoms / density);
        } else {
            double L = Math.pow(4.0 / density, 1.0 / 3.0);
            int n = (int) Math.round(Math.pow(numAtoms / 4, 1.0 / 3.0));
            boundaryHarmonic = new BoundaryRectangularPeriodic(space, n * L);
        }
        boxHarmonic.setBoundary(boundaryHarmonic);

        CoordinateDefinitionLeaf coordinateDefinitionHarmonic = new CoordinateDefinitionLeaf(boxHarmonic, primitive, basis, space);
        coordinateDefinitionHarmonic.initializeCoordinates(nCells);

        normalModes = new NormalModesFromFile(filename, space.D());
        normalModes.setHarmonicFudge(harmonicFudge);
        /*
         * nuke this line if it is DB
         */
        normalModes.setTemperature(temperature);

        WaveVectorFactory waveVectorFactory = normalModes.getWaveVectorFactory();
        waveVectorFactory.makeWaveVectors(boxHarmonic);
        move.setOmegaSquared(normalModes.getOmegaSquared());
        move.setEigenVectors(normalModes.getEigenvectors());
        move.setWaveVectors(waveVectorFactory.getWaveVectors());
        move.setWaveVectorCoefficients(waveVectorFactory.getCoefficients());
        move.setCoordinateDefinition(coordinateDefinitionHarmonic);
        move.setTemperature(temperature);

        move.setBox(boxHarmonic);

        // OVERLAP
        integratorOverlap = new IntegratorOverlap(new IntegratorBox[]{integratorHarmonic, integratorTarget});
        meterHarmonicEnergy = new MeterHarmonicEnergy(coordinateDefinitionTarget, normalModes);
        MeterBoltzmannTarget meterTarget = new MeterBoltzmannTarget(new MeterPotentialEnergyFromIntegrator(integratorTarget), meterHarmonicEnergy);
        meterTarget.setTemperature(temperature);
        meterTarget.setLatticeEnergy(latticeEnergy);
        meters[1] = meterTarget;
        setAccumulator(new AccumulatorVirialOverlapSingleAverage(10, 11, false), 1);

        MeterBoltzmannHarmonic meterHarmonic = new MeterBoltzmannHarmonic(move, potentialMasterTarget);
        meterHarmonic.setTemperature(temperature);
        meterHarmonic.setLatticeEnergy(latticeEnergy);
        meters[0] = meterHarmonic;
        setAccumulator(new AccumulatorVirialOverlapSingleAverage(10, 11, true), 0);

        setRefPref(1.0, 30);

        activityIntegrate = new ActivityIntegrate(integratorOverlap);

        getController().addAction(activityIntegrate);
    }

    /**
     * @param args filename containing simulation parameters
     * @see SimOverlapSoftSphereReweighting.SimOverlapParam
     */
    public static void main(String[] args) {

        //set up simulation parameters
        SimOverlapParam params = new SimOverlapParam();
        String inputFilename = null;
        if (args.length > 0) {
            inputFilename = args[0];
        }
        if (inputFilename != null) {
            ReadParameters readParameters = new ReadParameters(inputFilename, params);
            readParameters.readParameters();
        }
        double density = params.density/1000;
        int exponentN = params.exponentN;
        long numSteps = params.numSteps;
        final int numMolecules = params.numMolecules;
        double harmonicFudge = params.harmonicFudge;
        double temperature = params.temperature;
        int D = params.D;
        String filename = params.filename;
        if (filename.length() == 0) {
        	System.err.println("Need input files!!!");
            filename = "CB_FCC_n"+exponentN+"_T"+ (int)Math.round(temperature*10);
        }
        //String refFileName = args.length > 0 ? filename+"_ref" : null;
        String refFileName = filename+"_ref";


        System.out.println("Running "+(D==1 ? "1D" : (D==3 ? "FCC" : "2D hexagonal")) +" soft sphere overlap simulation");
        System.out.println(numMolecules+" atoms at density "+density+" and temperature "+temperature);
        System.out.println("exponent N: "+ exponentN +" and harmonic fudge: "+harmonicFudge);
        System.out.println((numSteps/1000)+" total steps of 1000");
        System.out.println("output data to "+filename);

        //instantiate simulation
        final SimOverlapSoftSphereReweighting sim = new SimOverlapSoftSphereReweighting(Space.getInstance(D), numMolecules, density, temperature, filename, harmonicFudge, exponentN);

        //start simulation
        sim.integratorOverlap.setNumSubSteps(1000);
        numSteps /= 1000;

        sim.initRefPref(refFileName, numSteps/20);
        if (Double.isNaN(sim.refPref) || sim.refPref == 0 || Double.isInfinite(sim.refPref)) {
            throw new RuntimeException("Simulation failed to find a valid ref pref");
        }
        System.out.flush();

        sim.equilibrate(refFileName, numSteps/10);
        if (Double.isNaN(sim.refPref) || sim.refPref == 0 || Double.isInfinite(sim.refPref)) {
            throw new RuntimeException("Simulation failed to find a valid ref pref");
        }

        System.out.println("equilibration finished");
        System.out.flush();

        // Harmonic to Target Sampling
        MeterWorkHarmonicTargetandReweighting meterWorkHarmonic =
                new MeterWorkHarmonicTargetandReweighting(sim.move, sim.potentialMasterTarget, sim.refPref);
        meterWorkHarmonic.setLatticeEnergy(sim.latticeEnergy);
        meterWorkHarmonic.setTemperature(temperature);

        DataFork dataForkHarmonic = new DataFork();
        DataPumpListener dataPumpHarmonic = new DataPumpListener(meterWorkHarmonic, dataForkHarmonic,1);

        AccumulatorAverageFixed dataAverageHarmonic = new AccumulatorAverageFixed(1);
        AccumulatorHistogram histogramHarmonicTarget = new AccumulatorHistogram(new HistogramSimple(4000, new DoubleRange(-200, 200)));
        dataForkHarmonic.addDataSink(dataAverageHarmonic);
        dataForkHarmonic.addDataSink(histogramHarmonicTarget);

        sim.integrators[0].getEventManager().addListener(dataPumpHarmonic);

        // Target to Harmonic Sampling
        MeterWorkTargetHarmonicandReweighting meterWorkTarget =
                new MeterWorkTargetHarmonicandReweighting(sim.integrators[1], sim.meterHarmonicEnergy, sim.refPref);
        meterWorkTarget.setLatticeEnergy(sim.latticeEnergy);

        DataFork dataForkTarget = new DataFork();
        DataPumpListener dataPumpTarget = new DataPumpListener(meterWorkTarget, dataForkTarget, 1000);

        AccumulatorAverageFixed dataAverageTarget = new AccumulatorAverageFixed(1);
        AccumulatorHistogram histogramTargetHarmonic = new AccumulatorHistogram(new HistogramSimple(4000, new DoubleRange(-200, 200)));
        dataForkTarget.addDataSink(dataAverageTarget);
        dataForkTarget.addDataSink(histogramTargetHarmonic);

        sim.integrators[1].getEventManager().addListener(dataPumpTarget);

        int totalCells = 1;
        for (int i=0; i<D; i++) {
            totalCells *= sim.nCells[i];
        }
        int basisSize = sim.basis.getScaledCoordinates().length;

        double  AHarmonic = CalcHarmonicA.doit(sim.normalModes, D, temperature, basisSize*totalCells);
        System.out.println("Harmonic-reference free energy, A: "+AHarmonic + " " + AHarmonic/numMolecules);

        sim.activityIntegrate.setMaxSteps(numSteps);
        sim.getController().actionPerformed();

        System.out.println("final reference optimal step frequency "+sim.integratorOverlap.getIdealRefStepFraction()
        		+" (actual: "+sim.integratorOverlap.getRefStepFraction()+")");
        double[] ratioAndError = sim.dsvo.getOverlapAverageAndError();
        double ratio = ratioAndError[0];
        double error = ratioAndError[1];
        double betaFAB = -Math.log(ratio);

        System.out.println("\nratio average: "+ratio+" ,error: "+error);
        System.out.println("free energy difference: "+(-temperature*Math.log(ratio))+" ,error: "+temperature*(error/ratio));
        System.out.println("target free energy: "+temperature*(AHarmonic-Math.log(ratio)));
        System.out.println("target free energy per particle: "+temperature*(AHarmonic-Math.log(ratio))/numMolecules);

        DataGroup allYourBase = (DataGroup)sim.accumulators[0].getData(sim.dsvo.minDiffLocation());
        double betaFAW = -Math.log(((DataDoubleArray) allYourBase.getData(AccumulatorAverage.AVERAGE.index)).getData()[1]);
        System.out.println("harmonic ratio average: " + ((DataDoubleArray) allYourBase.getData(AccumulatorAverage.AVERAGE.index)).getData()[1]
                + " stdev: " + ((DataDoubleArray) allYourBase.getData(AccumulatorAverage.STANDARD_DEVIATION.index)).getData()[1]
                + " error: " + ((DataDoubleArray) allYourBase.getData(AccumulatorAverage.ERROR.index)).getData()[1]);

        allYourBase = (DataGroup)sim.accumulators[1].getData(sim.dsvo.minDiffLocation());
        double betaFBW = -Math.log(((DataDoubleArray) allYourBase.getData(AccumulatorAverage.AVERAGE.index)).getData()[1]);
        System.out.println("target ratio average: " + ((DataDoubleArray) allYourBase.getData(AccumulatorAverage.AVERAGE.index)).getData()[1]
                + " stdev: " + ((DataDoubleArray) allYourBase.getData(AccumulatorAverage.STANDARD_DEVIATION.index)).getData()[1]
                + " error: " + ((DataDoubleArray) allYourBase.getData(AccumulatorAverage.ERROR.index)).getData()[1]);

        /*
         * Refer Wu & Kofke JCp 123,054103(2003) Eq (6)
         */
        double betaUAB = dataAverageHarmonic.getData().getValue(AccumulatorAverage.AVERAGE.index);
        double betaUBA = dataAverageTarget.getData().getValue(AccumulatorAverage.AVERAGE.index);

        double SAB = betaUAB - betaFAB;
        double SBA = betaUBA + betaFAB;

        double betaUAWf = meterWorkHarmonic.getBetaUAWf();  // < beta*U_AW>A
        double betaUAWr = meterWorkHarmonic.getBetaUAWr();  // < beta*U_AW>W

        double betaUBWf = meterWorkTarget.getBetaUBWf();	// < beta*U_BW>B
        double betaUBWr = meterWorkTarget.getBetaUBWr();	// < beta*U_BW>W

        double SAW =  betaUAWf - betaFAW;
        double SWA = - betaUAWr + betaFAW;

        double SBW = betaUBWf - betaFBW;
        double SWB = - betaUBWr + betaFBW;

        System.out.println("");
        System.out.println("SAB: "+ SAB + " ; betaUAB: " + betaUAB + " ;betaFAB: " + betaFAB);
        System.out.println("SBA: "+ SBA + " ; betaUBA: " + betaUBA + " ;betaFAB: " + betaFAB);

        System.out.println("");
        System.out.println("SAW: "+ SAW + " ; betaUAWf: " + betaUAWf + " ;betaFAW: " + betaFAW);
        System.out.println("SWA: "+ SWA + " ; betaUAWr: " + betaUAWr + " ;betaFAW: " + betaFAW);

        System.out.println("");
        System.out.println("SBW: "+ SBW + " ; betaUBWf: " + betaUBWf + " ;betaFBW: " + betaFBW);
        System.out.println("SWB: "+ SWB + " ; betaUBWr: " + betaUBWr + " ;betaFBW: " + betaFBW);

        /*
         * Harmonic -----> Target
         */
        // AB
        DataLogger dataLogger = new DataLogger();
        DataTableWriter dataTableWriter = new DataTableWriter();
        dataLogger.setFileName(filename + "_hist_HarmTarg");
        dataTableWriter.setIncludeHeader(false);

        dataLogger.setDataSink(dataTableWriter);
        dataLogger.putDataInfo(histogramHarmonicTarget.getDataInfo());
        dataLogger.putData(histogramHarmonicTarget.getData());
        dataLogger.closeFile();

        // AW
        dataLogger = new DataLogger();
        dataTableWriter = new DataTableWriter();
        dataLogger.setFileName(filename + "_hist_HarmBenn");
        dataTableWriter.setIncludeHeader(false);

        dataLogger.setDataSink(dataTableWriter);
        dataLogger.putDataInfo(meterWorkHarmonic.getDataInfoHistogramBetaUAWf());
        dataLogger.putData(meterWorkHarmonic.getDataHistogramBetaUAWf());
        dataLogger.closeFile();

        // WA
        dataLogger = new DataLogger();
        dataTableWriter = new DataTableWriter();
        dataLogger.setFileName(filename + "_hist_BennHarm");
        dataTableWriter.setIncludeHeader(false);

        dataLogger.setDataSink(dataTableWriter);
        dataLogger.putDataInfo(meterWorkHarmonic.getDataInfoHistogramBetaUAWr());
        dataLogger.putData(meterWorkHarmonic.getDataHistogramBetaUAWr());
        dataLogger.closeFile();

        /*
         * Target -----> Harmonic
         */
        //BA
        dataLogger = new DataLogger();
        dataTableWriter = new DataTableWriter();
        dataLogger.setFileName(filename + "_hist_TargHarm");
        dataTableWriter.setIncludeHeader(false);

        dataLogger.setDataSink(dataTableWriter);
        dataLogger.putDataInfo(histogramTargetHarmonic.getDataInfo());
        dataLogger.putData(histogramTargetHarmonic.getData());
        dataLogger.closeFile();

        //BW
        dataLogger = new DataLogger();
        dataTableWriter = new DataTableWriter();
        dataLogger.setFileName(filename + "_hist_TargBenn");
        dataTableWriter.setIncludeHeader(false);

        dataLogger.setDataSink(dataTableWriter);
        dataLogger.putDataInfo(meterWorkTarget.getDataInfoHistogramBetaUBWf());
        dataLogger.putData(meterWorkTarget.getDataHistogramBetaUBWf());
        dataLogger.closeFile();

        // WB
        dataLogger = new DataLogger();
        dataTableWriter = new DataTableWriter();
        dataLogger.setFileName(filename + "_hist_BennTarg");
        dataTableWriter.setIncludeHeader(false);

        dataLogger.setDataSink(dataTableWriter);
        dataLogger.putDataInfo(meterWorkTarget.getDataInfoHistogramBetaUBWr());
        dataLogger.putData(meterWorkTarget.getDataHistogramBetaUBWr());
        dataLogger.closeFile();

    }

    public void setRefPref(double refPrefCenter, double span) {
        refPref = refPrefCenter;
        accumulators[0].setBennetParam(refPrefCenter, span);
        accumulators[1].setBennetParam(refPrefCenter, span);

    }

    public void setAccumulator(AccumulatorVirialOverlapSingleAverage newAccumulator, int iBox) {

        accumulators[iBox] = newAccumulator;

        newAccumulator.setBlockSize(200);

        if (accumulatorPumps[iBox] == null) {
            accumulatorPumps[iBox] = new DataPumpListener(meters[iBox], newAccumulator);
            integrators[iBox].getEventManager().addListener(accumulatorPumps[iBox]);
            if (iBox == 1) {
                if (boxTarget.getMoleculeList().size() == 32) {
                    accumulatorPumps[iBox].setInterval(100);

                } else if (boxTarget.getMoleculeList().size() == 108) {

                    accumulatorPumps[iBox].setInterval(300);
                } else
                    accumulatorPumps[iBox].setInterval(boxTarget.getMoleculeList().size());
            }
        } else {
            accumulatorPumps[iBox].setDataSink(newAccumulator);
        }
        if (integratorOverlap != null && accumulators[0] != null && accumulators[1] != null) {
            dsvo = new DataSourceVirialOverlap(accumulators[0], accumulators[1]);
            integratorOverlap.setReferenceFracSource(dsvo);
        }
    }

    public void setRefPref(double newRefPref) {
        System.out.println("setting ref pref (explicitly) to " + newRefPref);
        setAccumulator(new AccumulatorVirialOverlapSingleAverage(1, true), 0);
        setAccumulator(new AccumulatorVirialOverlapSingleAverage(1, false), 1);
        setRefPref(newRefPref, 1);
    }

    public void initRefPref(String fileName, long initSteps) {
        // refPref = -1 indicates we are searching for an appropriate value
        refPref = -1.0;
        if (fileName != null) {
            try {
                FileReader fileReader = new FileReader(fileName);
                BufferedReader bufReader = new BufferedReader(fileReader);
                String refPrefString = bufReader.readLine();
                refPref = Double.parseDouble(refPrefString);
                bufReader.close();
                fileReader.close();
                System.out.println("setting ref pref (from file) to " + refPref);
                setAccumulator(new AccumulatorVirialOverlapSingleAverage(1, true), 0);
                setAccumulator(new AccumulatorVirialOverlapSingleAverage(1, false), 1);
                setRefPref(refPref, 1);
            } catch (IOException e) {
                // file not there, which is ok.
            }
        }

        if (refPref == -1) {
            // equilibrate off the lattice to avoid anomolous contributions
            activityIntegrate.setMaxSteps(initSteps / 2);
            getController().actionPerformed();
            getController().reset();
            System.out.println("target equilibration finished");

            setAccumulator(new AccumulatorVirialOverlapSingleAverage(41, true), 0);
            setAccumulator(new AccumulatorVirialOverlapSingleAverage(41, false), 1);
            setRefPref(1, 200);
            activityIntegrate.setMaxSteps(initSteps);
            getController().actionPerformed();
            getController().reset();

            int newMinDiffLoc = dsvo.minDiffLocation();
            refPref = accumulators[0].getBennetAverage(newMinDiffLoc)
                    / accumulators[1].getBennetAverage(newMinDiffLoc);
            if (Double.isNaN(refPref) || refPref == 0 || Double.isInfinite(refPref)) {
                throw new RuntimeException("Simulation failed to fiq" +
                        "nd a valid ref pref");
            }
            System.out.println("setting ref pref to " + refPref);

            setAccumulator(new AccumulatorVirialOverlapSingleAverage(11, true), 0);
            setAccumulator(new AccumulatorVirialOverlapSingleAverage(11, false), 1);
            setRefPref(refPref, 5);

            // set refPref back to -1 so that later on we know that we've been looking for
            // the appropriate value
            refPref = -1;
            getController().reset();
        }

    }

    public void equilibrate(String fileName, long initSteps) {
        // run a short simulation to get reasonable MC Move step sizes and
        // (if needed) narrow in on a reference preference
        activityIntegrate.setMaxSteps(initSteps);

        for (int i = 0; i < 2; i++) {
            if (integrators[i] instanceof IntegratorMC)
                ((IntegratorMC) integrators[i]).getMoveManager().setEquilibrating(true);
        }
        getController().actionPerformed();
        getController().reset();
        for (int i = 0; i < 2; i++) {
            if (integrators[i] instanceof IntegratorMC)
                ((IntegratorMC) integrators[i]).getMoveManager().setEquilibrating(false);
        }

        if (refPref == -1) {
            int newMinDiffLoc = dsvo.minDiffLocation();
            refPref = accumulators[0].getBennetAverage(newMinDiffLoc)
                    / accumulators[1].getBennetAverage(newMinDiffLoc);
            System.out.println("setting ref pref to " + refPref + " (" + newMinDiffLoc + ")");
            setAccumulator(new AccumulatorVirialOverlapSingleAverage(1, true), 0);

            System.out.println("block size (equilibrate) " + accumulators[0].getBlockSize());

            setAccumulator(new AccumulatorVirialOverlapSingleAverage(1, false), 1);
            setRefPref(refPref, 1);
            if (fileName != null) {
                try {
                    FileWriter fileWriter = new FileWriter(fileName);
                    BufferedWriter bufWriter = new BufferedWriter(fileWriter);
                    bufWriter.write(String.valueOf(refPref) + "\n");
                    bufWriter.close();
                    fileWriter.close();
                } catch (IOException e) {
                    throw new RuntimeException("couldn't write to refpref file");
                }
            }
        } else {
            dsvo.reset();
        }
    }

    /**
     * Inner class for parameters understood by the HSMD3D constructor
     */
    public static class SimOverlapParam extends ParameterBase {
        public int numMolecules =32;
        public double density = 1256;
        public int exponentN = 12;
        public int D = 3;
        public long numSteps = 100000;
        public double harmonicFudge = 1;
        public String filename = "CB_FCC_n12_T01";
        public double temperature = 0.1;
    }
}
