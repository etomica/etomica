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
import etomica.lattice.crystal.Basis;
import etomica.lattice.crystal.BasisCubicFcc;
import etomica.lattice.crystal.Primitive;
import etomica.lattice.crystal.PrimitiveCubic;
import etomica.integrator.IntegratorListenerAction;
import etomica.math.DoubleRange;
import etomica.nbr.list.PotentialMasterList;
import etomica.overlap.IntegratorOverlap;
import etomica.potential.P2SoftSphere;
import etomica.potential.P2SoftSphericalTruncated;
import etomica.potential.Potential2SoftSpherical;
import etomica.simulation.Simulation;
import etomica.space.Boundary;
import etomica.space.BoundaryRectangularPeriodic;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.space3d.Vector3D;
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
 * The original Bennett's Overlapping Sampling Simulation
 *  - used to check for the computation time
 * 
 * @author Nancy Cribbin of TBT's SimOverlapSoftSphere
 */


public class SimOverlapSSnxy extends Simulation {

    public IntegratorOverlap integratorOverlap;
    public DataSourceVirialOverlap dsvo;
    public IntegratorBox[] integrators;
    public ActivityIntegrate activityIntegrate;
    public Box boxTarget, boxHarmonic;
    public Boundary boundaryTarget, boundaryHarmonic;
    public int[] nCells;
    public Basis basis;
    public NormalModes normalModes;
    public Primitive primitive, primitiveUnitCell;
    public double refPref;
    public AccumulatorVirialOverlapSingleAverage[] accumulators;
    public DataPump[] accumulatorPumps;
    public IDataSource[] meters;
    public String fname;
    protected MCMoveHarmonic move;
    protected MCMoveAtomCoupled atomMove;
    protected PotentialMasterList potentialMasterTarget;
    protected MeterHarmonicEnergy meterHarmonicEnergy;
    protected double latticeEnergy;
    
    public SimOverlapSSnxy(Space _space, int numAtoms, double density, double 
            temperature, String filename, double harmonicFudge, int exponent, 
            int[] shape) {
        super(_space);
        this.fname = filename;
        
        integrators = new IntegratorBox[2];
        accumulatorPumps = new DataPump[2];
        meters = new IDataSource[2];
        accumulators = new AccumulatorVirialOverlapSingleAverage[2];

        SpeciesSpheresMono species = new SpeciesSpheresMono(this, space);
        addSpecies(species);

        // TARGET
        boxTarget = new Box(space);
        addBox(boxTarget);
        boxTarget.setNMolecules(species, numAtoms);
        
        double primitiveLength = Math.pow(4.0 / density, 1.0 / 3.0);
        primitive = new PrimitiveCubic(space);
        basis = new BasisCubicFcc();
        primitive.scaleSize(primitiveLength);
        
        boundaryTarget = new BoundaryRectangularPeriodic(space, 1.0);
        Vector edges = new Vector3D();
        double[] lengths = new double[3];
        lengths[0] = shape[0]*primitiveLength;
        lengths[1] = shape[1]*primitiveLength;
        lengths[2] = shape[2]*primitiveLength;
        edges.E(lengths);
        boundaryTarget.setBoxSize(edges);
        boxTarget.setBoundary(boundaryTarget);
        
        potentialMasterTarget = new PotentialMasterList(this, space);
        double neighborRange = 0.495 * lengths[0];
        Potential2SoftSpherical potentialBase = new P2SoftSphere(space, 1.0, 
                1.0, exponent);
        P2SoftSphericalTruncated potential = new P2SoftSphericalTruncated(
                space, potentialBase, neighborRange);
        potentialMasterTarget.addPotential(potential, new AtomType[]{
                species.getLeafType(), species.getLeafType()});
        potentialMasterTarget.lrcMaster().setEnabled(false);
        potentialMasterTarget.setRange(neighborRange);
        potentialMasterTarget.getNeighborManager(boxTarget).reset();
        
        CoordinateDefinitionLeaf coordinateDefinitionTarget = new 
                CoordinateDefinitionLeaf(boxTarget, primitive, basis, space);
        coordinateDefinitionTarget.initializeCoordinates(shape);
        coordinateDefinitionTarget = new CoordinateDefinitionLeaf(boxTarget, 
                primitive, basis, space);
        coordinateDefinitionTarget.initializeCoordinates(shape);

        IntegratorMC integratorTarget = new IntegratorMC(potentialMasterTarget, 
                getRandom(), temperature);
        integratorTarget.setBox(boxTarget);
        integrators[1] = integratorTarget;
        
        atomMove = new MCMoveAtomCoupled(potentialMasterTarget, new 
                MeterPotentialEnergy(potentialMasterTarget), getRandom(), space);
        atomMove.setStepSize(0.1);
        atomMove.setStepSizeMax(0.5);
        atomMove.setDoExcludeNonNeighbors(true);
        atomMove.setPotential(potential);
        integratorTarget.getMoveManager().addMCMove(atomMove);
        ((MCMoveStepTracker)atomMove.getTracker()).setNoisyAdjustment(true);
        
        MeterPotentialEnergy meterPE = new MeterPotentialEnergy(potentialMasterTarget);
        meterPE.setBox(boxTarget);
        latticeEnergy = meterPE.getDataAsScalar();
        System.out.println("lattice energy of SS system: " + latticeEnergy);
       
        // HARMONIC
        boxHarmonic = new Box(space);
        addBox(boxHarmonic);
        boxHarmonic.setNMolecules(species, numAtoms);

        IntegratorMC integratorHarmonic = new IntegratorMC(null, random, 1.0); 
        move = new MCMoveHarmonic(getRandom());
        integratorHarmonic.getMoveManager().addMCMove(move);
        integrators[0] = integratorHarmonic;
        
        boundaryHarmonic = new BoundaryRectangularPeriodic(space, 1.0);
        boundaryHarmonic.setBoxSize(edges);
        boxHarmonic.setBoundary(boundaryHarmonic);
        CoordinateDefinitionLeaf coordinateDefinitionHarmonic = new 
                CoordinateDefinitionLeaf(boxHarmonic, primitive, basis, space);
        coordinateDefinitionHarmonic.initializeCoordinates(shape);
        
        String inFile = "inputSSDB_"+numAtoms;
        normalModes = new NormalModesFromFile(inFile, space.D());
        normalModes.setHarmonicFudge(harmonicFudge);
        /*
         * nuke this line if it is DB
         */
        //normalModes.setTemperature(temperature);
        
        WaveVectorFactory waveVectorFactory = normalModes.getWaveVectorFactory();
        waveVectorFactory.makeWaveVectors(boxHarmonic);
        move.setOmegaSquared(normalModes.getOmegaSquared());
        move.setEigenVectors(normalModes.getEigenvectors());
        move.setWaveVectors(waveVectorFactory.getWaveVectors());
        move.setWaveVectorCoefficients(waveVectorFactory.getCoefficients());
        move.setCoordinateDefinition(coordinateDefinitionHarmonic);
        move.setTemperature(temperature);
        
        move.setBox(boxHarmonic);
        
        integratorHarmonic.setBox(boxHarmonic);
        
        if (potentialMasterTarget instanceof PotentialMasterList) {
            // find neighbors now.  Don't hook up NeighborListManager (neighbors won't change)
            potentialMasterTarget.
                    getNeighborManager(boxHarmonic).reset();
        }
        
        // OVERLAP
        integratorOverlap = new IntegratorOverlap(new IntegratorBox[]{
                integratorHarmonic, integratorTarget});
        meterHarmonicEnergy = new MeterHarmonicEnergy(coordinateDefinitionTarget,
                normalModes);
        MeterBoltzmannTarget meterTarget = new MeterBoltzmannTarget(new 
                MeterPotentialEnergyFromIntegrator(integratorTarget), meterHarmonicEnergy);
        meterTarget.setFrac(1.0);
        meterTarget.setTemperature(temperature);
        meterTarget.setLatticeEnergy(latticeEnergy);
        meters[1] = meterTarget;
        setAccumulator(new AccumulatorVirialOverlapSingleAverage(10, 11, false), 1);
        
        MeterBoltzmannHarmonic meterHarmonic = new MeterBoltzmannHarmonic(move, 
                potentialMasterTarget);
        meterHarmonic.setFrac(1.0);
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
     * @see SimOverlapSoftSphere.SimOverlapParam
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
        double density = params.density/10000;
        int exponentN = params.exponentN;
        long numSteps = params.numSteps;
        final int numMolecules = params.numMolecules;
        double harmonicFudge = params.harmonicFudge;
        double temperature = params.temperature;
        int D = params.D;
        int[] shape = params.shape;
        String filename = params.filename;
        if (filename.length() == 0) {
            System.err.println("Need input files!!!");
            filename = "CB_FCC_n"+exponentN+"_T"+ (int)Math.round(temperature*10);
        }
        //String refFileName = args.length > 0 ? filename+"_ref" : null;
        String refFileName = filename+"_ref";


        System.out.println("Running "+(D==1 ? "1D" : (D==3 ? "FCC" : "2D " +
                "hexagonal")) +" soft sphere overlap simulation");
        System.out.println(numMolecules+" atoms at density "+density+" and " +
                "temperature "+temperature);
        System.out.println("exponent N: " + exponentN + " and harmonic fudge: " +
                harmonicFudge);
        System.out.println((numSteps/1000)+" total steps of 1000");
        System.out.println("output data to "+filename);

        //instantiate simulation
        final SimOverlapSSnxy sim = new SimOverlapSSnxy(Space.getInstance(D),
                numMolecules, density, temperature, filename, harmonicFudge,
                exponentN, shape);

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


        /*
         * To quantify the relative entropy
         */

        // Harmonic to Bennett Sampling
        MeterWorkHarmonicBennet meterWorkHarmonic =
                new MeterWorkHarmonicBennet(sim.move, sim.potentialMasterTarget,
                sim.refPref);
        meterWorkHarmonic.setLatticeEnergy(sim.latticeEnergy);
        meterWorkHarmonic.setTemperature(temperature);

        DataFork dataForkHarmonic = new DataFork();
        DataPumpListener dataPumpHarmonic = new DataPumpListener(
                meterWorkHarmonic, dataForkHarmonic,1);

        AccumulatorAverageFixed dataAverageHarmonic = new AccumulatorAverageFixed(1);
        AccumulatorHistogram histogramHarmonicTarget = new AccumulatorHistogram(
                new HistogramSimple(4000, new DoubleRange(-200, 200)));
        dataForkHarmonic.addDataSink(dataAverageHarmonic);
        dataForkHarmonic.addDataSink(histogramHarmonicTarget);
        sim.integrators[0].getEventManager().addListener(dataPumpHarmonic);

        // Target to Bennett Sampling
        MeterWorkTargetBennet meterWorkTarget = new MeterWorkTargetBennet(
                sim.integrators[1], sim.meterHarmonicEnergy, sim.refPref);
        meterWorkTarget.setLatticeEnergy(sim.latticeEnergy);

        DataFork dataForkTarget = new DataFork();
        DataPumpListener dataPumpTarget = new DataPumpListener(meterWorkTarget,
                dataForkTarget, 1000);

        AccumulatorAverageFixed dataAverageTarget = new AccumulatorAverageFixed(1);
        AccumulatorHistogram histogramTargetHarmonic = new AccumulatorHistogram(
                new HistogramSimple(4000, new DoubleRange(-200, 200)));
        dataForkTarget.addDataSink(dataAverageTarget);
        dataForkTarget.addDataSink(histogramTargetHarmonic);
        sim.integrators[1].getEventManager().addListener(dataPumpTarget);


        final long startTime = System.currentTimeMillis();
        System.out.println("Start Time: " + startTime);

        sim.activityIntegrate.setMaxSteps(numSteps);
        sim.getController().actionPerformed();

        int totalCells = 1;
        for (int i=0; i<D; i++) {
            totalCells *= shape[i];
        }
        double AHarmonic = CalcHarmonicA.doit(sim.normalModes, D, temperature,
                numMolecules);
        System.out.println("Harmonic-reference free energy, A: " + AHarmonic +
                " " + AHarmonic/numMolecules);
        System.out.println(" ");

        System.out.println("final reference optimal step frequency " +
                sim.integratorOverlap.getIdealRefStepFraction() + " (actual: " +
                sim.integratorOverlap.getRefStepFraction() + ")");

        double[] ratioAndError = sim.dsvo.getOverlapAverageAndError();
        double ratio = ratioAndError[0];
        double error = ratioAndError[1];

        System.out.println("\nratio average: " + ratio + " error: " + error);
        System.out.println("free energy difference: " + (-temperature*Math.log(ratio))
                +" error: "+temperature*(error/ratio));
        System.out.println("target free energy: "+(AHarmonic-temperature*
                Math.log(ratio)));
        System.out.println("target free energy per particle: "+ (AHarmonic-
                temperature * Math.log(ratio)) / numMolecules
                +" error: "+temperature*(error/ratio)/numMolecules);
        DataGroup allYourBase = (DataGroup)sim.accumulators[0].getData(
                sim.dsvo.minDiffLocation());
        System.out.println("harmonic ratio average: "+((DataDoubleArray)
                allYourBase.getData(AccumulatorAverage.AVERAGE.index)).getData()[1]
                +" stdev: "+((DataDoubleArray)allYourBase.getData(
                AccumulatorAverage.STANDARD_DEVIATION.index)).getData()[1]
                +" error: "+((DataDoubleArray)allYourBase.getData(
                AccumulatorAverage.ERROR.index)).getData()[1]);

        allYourBase = (DataGroup)sim.accumulators[1].getData(sim.dsvo.minDiffLocation());
        System.out.println("target ratio average: "+((DataDoubleArray)
                allYourBase.getData(AccumulatorAverage.AVERAGE.index)).getData()[1]
                +" stdev: "+((DataDoubleArray)allYourBase.getData(
                AccumulatorAverage.STANDARD_DEVIATION.index)).getData()[1]
                +" error: "+((DataDoubleArray)allYourBase.getData(
                AccumulatorAverage.ERROR.index)).getData()[1]);

        long endTime = System.currentTimeMillis();
        System.out.println("End Time: " + endTime);
        System.out.println("Time taken: " + (endTime - startTime));
    }

    public void setRefPref(double refPrefCenter, double span) {
        refPref = refPrefCenter;
        accumulators[0].setBennetParam(refPrefCenter, span);
        accumulators[1].setBennetParam(refPrefCenter, span);
    }

    public void setAccumulator(AccumulatorVirialOverlapSingleAverage
                                       newAccumulator, int iBox) {

        accumulators[iBox] = newAccumulator;

        newAccumulator.setBlockSize(200); // setting the block size = 300

        if (accumulatorPumps[iBox] == null) {
            accumulatorPumps[iBox] = new DataPump(meters[iBox], newAccumulator);
            IntegratorListenerAction pumpListener = new
                    IntegratorListenerAction(accumulatorPumps[iBox]);
            integrators[iBox].getEventManager().addListener(pumpListener);
            if (iBox == 1) {
                if (boxTarget.getMoleculeList().getMoleculeCount() == 32) {
                    pumpListener.setInterval(100);
                } else if (boxTarget.getMoleculeList().getMoleculeCount() == 108) {
                    pumpListener.setInterval(300);
                } else
                    pumpListener.setInterval(boxTarget.getMoleculeList().
                            getMoleculeCount());
            }
        } else {
            accumulatorPumps[iBox].setDataSink(newAccumulator);
        }
        if (integratorOverlap != null && accumulators[0] != null &&
                accumulators[1] != null) {
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
            if (integrators[i] instanceof IntegratorMC) ((IntegratorMC)
                    integrators[i]).getMoveManager().setEquilibrating(true);
        }
        getController().actionPerformed();
        getController().reset();
        for (int i = 0; i < 2; i++) {
            if (integrators[i] instanceof IntegratorMC) ((IntegratorMC)
                    integrators[i]).getMoveManager().setEquilibrating(false);
        }

        if (refPref == -1) {
            int newMinDiffLoc = dsvo.minDiffLocation();
            refPref = accumulators[0].getBennetAverage(newMinDiffLoc)
                    / accumulators[1].getBennetAverage(newMinDiffLoc);
            System.out.println("setting ref pref to " + refPref + " (" +
                    newMinDiffLoc + ")");
            setAccumulator(new AccumulatorVirialOverlapSingleAverage(1, true), 0);

            System.out.println("block size (equilibrate) " +
                    accumulators[0].getBlockSize());

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
        public int numMolecules = 32;
        public double density = 11964;
        public int exponentN = 12;
        public int D = 3;
        public long numSteps = 1000000;
        public double harmonicFudge = 1;
        public String filename = "inputSSDB32T01";
        public double temperature = 0.01;
        public int[] shape = {2, 2, 2};
    }
}
