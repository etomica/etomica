/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.models.oneDHardRods;

import etomica.action.activity.ActivityIntegrate;
import etomica.atom.AtomType;
import etomica.box.Box;
import etomica.data.AccumulatorRatioAverageCovariance;
import etomica.data.DataPump;
import etomica.data.IDataSource;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.data.types.DataDoubleArray;
import etomica.data.types.DataGroup;
import etomica.exception.ConfigurationOverlapException;
import etomica.integrator.IntegratorMC;
import etomica.lattice.crystal.BasisMonatomic;
import etomica.lattice.crystal.Primitive;
import etomica.lattice.crystal.PrimitiveCubic;
import etomica.integrator.IntegratorListenerAction;
import etomica.math.SpecialFunctions;
import etomica.nbr.list.PotentialMasterList;
import etomica.normalmode.*;
import etomica.overlap.IntegratorOverlap;
import etomica.potential.P2HardSphere;
import etomica.potential.Potential2;
import etomica.potential.Potential2HardSpherical;
import etomica.simulation.Simulation;
import etomica.space.Boundary;
import etomica.space.BoundaryRectangularPeriodic;
import etomica.space.Space;
import etomica.species.SpeciesSpheresMono;
import etomica.units.dimensions.Null;
import etomica.util.ParameterBase;
import etomica.util.ReadParameters;
import etomica.virial.overlap.AccumulatorVirialOverlapSingleAverage;
import etomica.virial.overlap.DataSourceVirialOverlap;

import java.io.*;

public class SimOverlapMultipleWaveVectors extends Simulation {

    private static final String APP_NAME = "SimOverlapMultipleWaveVectors";
    public IntegratorOverlap integratorSim; //integrator for the whole simulation
    public DataSourceVirialOverlap dsvo;
    public BasisMonatomic basis;
    public AccumulatorVirialOverlapSingleAverage[] accumulators;
    public DataPump[] accumulatorPumps;
    public IDataSource[] meters;
    public Box boxTarget, boxRef;
    public Boundary boundaryTarget, boundaryRef;
    Primitive primitive;
    int[] nCells;
    NormalModes1DHR nm;
    double bennettParam;       //adjustable parameter - Bennett's parameter
    ActivityIntegrate activityIntegrate;
    IntegratorMC[] integrators;
    MCMoveChangeMultipleWV changeMove;
    MCMoveCompareMultipleWV compareMove;
    MeterPotentialEnergy meterAinB, meterAinA;
    MeterCompareMultipleWVBrute meterBinA, meterBinB;
    
    
    public SimOverlapMultipleWaveVectors(Space _space, int numAtoms, double 
            density, double temperature, String filename, double harmonicFudge,
            int[] compWV, int[] chbleWV) {
        super(_space);

        System.out.println("Running " + SimOverlapMultipleWaveVectors.APP_NAME);

//        long seed = 2;
//        System.out.println("Seed explicitly set to " + seed);
//        IRandom rand = new RandomNumberGenerator(seed);
//        this.setRandom(rand);

        //Set up some of the joint stuff
        SpeciesSpheresMono species = new SpeciesSpheresMono(this, space);
        addSpecies(species);

        integrators = new IntegratorMC[2];
        accumulatorPumps = new DataPump[2];
        meters = new IDataSource[2];
        accumulators = new AccumulatorVirialOverlapSingleAverage[2];

        basis = new BasisMonatomic(space);


//TARGET    
        // Set up target system - A, 1, hard rod
        PotentialMasterList potentialMasterTarget = new PotentialMasterList(
                this, space);
        boxTarget = this.makeBox();
        boxTarget.setNMolecules(species, numAtoms);

        Potential2 p2 = new P2HardSphere(space, 1.0, true);
        p2 = new P2XOrder(space, (Potential2HardSpherical) p2);
        p2.setBox(boxTarget);
        potentialMasterTarget.addPotential(p2, new AtomType[]{
                species.getLeafType(), species.getLeafType()});

        primitive = new PrimitiveCubic(space, 1.0 / density);
        boundaryTarget = new BoundaryRectangularPeriodic(space, numAtoms / density);
        nCells = new int[]{numAtoms};
        boxTarget.setBoundary(boundaryTarget);

        CoordinateDefinitionLeaf coordinateDefinitionTarget = new
                CoordinateDefinitionLeaf(boxTarget, primitive, space);
        coordinateDefinitionTarget.initializeCoordinates(nCells);

        double neighborRange = 1.01 / density;
        potentialMasterTarget.setRange(neighborRange);
        // Find neighbors now.  Don't hook up the NieghborListManager since the
        //  neighbors won't change.
        potentialMasterTarget.getNeighborManager(boxTarget).reset();

        IntegratorMC integratorTarget = new IntegratorMC(potentialMasterTarget,
                random, temperature, boxTarget);
        integrators[1] = integratorTarget;

        nm = new NormalModes1DHR(boundaryTarget, numAtoms);
        nm.setHarmonicFudge(harmonicFudge);
        nm.setTemperature(temperature);

        WaveVectorFactory waveVectorFactoryTarget = nm.getWaveVectorFactory();
        waveVectorFactoryTarget.makeWaveVectors(boxTarget);
        int wvflength = waveVectorFactoryTarget.getWaveVectors().length;
        System.out.println("We have " + wvflength + " wave vectors.");
        System.out.println("Wave Vector Coefficients:");
        for (int i = 0; i < wvflength; i++) {
            System.out.println(i + " " + waveVectorFactoryTarget.getCoefficients()[i]);
        }

        changeMove = new MCMoveChangeMultipleWV(potentialMasterTarget, random);
        integratorTarget.getMoveManager().addMCMove(changeMove);
        changeMove.setWaveVectors(waveVectorFactoryTarget.getWaveVectors());
        changeMove.setWaveVectorCoefficients(waveVectorFactoryTarget.getCoefficients());
        changeMove.setEigenVectors(nm.getEigenvectors());
        changeMove.setCoordinateDefinition(coordinateDefinitionTarget);
        changeMove.setBox(boxTarget);
        changeMove.setStepSizeMin(0.001);
        changeMove.setStepSize(0.01);
        changeMove.setOmegaSquared(nm.getOmegaSquared());

        meterAinA = new MeterPotentialEnergy(potentialMasterTarget, boxTarget);

        meterBinA = new MeterCompareMultipleWVBrute("meterBinA",
                potentialMasterTarget, coordinateDefinitionTarget, boxTarget);
        meterBinA.setEigenVectors(nm.getEigenvectors());
        meterBinA.setOmegaSquared(nm.getOmegaSquared());
        meterBinA.setTemperature(temperature);
        meterBinA.setWaveVectorCoefficients(waveVectorFactoryTarget.getCoefficients());
        meterBinA.setWaveVectors(waveVectorFactoryTarget.getWaveVectors());

        MeterOverlap meterOverlapInA = new MeterOverlap("meterOverlapInA", Null.DIMENSION,
                meterAinA, meterBinA, temperature);
        meters[1] = meterOverlapInA;

        potentialMasterTarget.getNeighborManager(boxTarget).reset();

//        meterBinA.getSingle().setCoordinateDefinition(coordinateDefinitionTarget);
//        meterBinA.getSingle().setEigenVectors(nm.getEigenvectors(boxTarget));
//        meterBinA.getSingle().setOmegaSquared(nm.getOmegaSquared(boxTarget));
//        meterBinA.getSingle().setTemperature(temperature);
//        meterBinA.getSingle().setWaveVectorCoefficients(waveVectorFactoryTarget.getCoefficients());
//        meterBinA.getSingle().setWaveVectors(waveVectorFactoryTarget.getWaveVectors());
//        meterBinA.setA(true);

//        singleBinA = new MeterCompareSingleModeBrute("singleBinA", potentialMasterTarget, coordinateDefinitionTarget, boxTarget);
//        singleBinA.setEigenVectors(nm.getEigenvectors(boxTarget));
//        singleBinA.setOmegaSquared(nm.getOmegaSquared(boxTarget));
//        singleBinA.setTemperature(temperature);
//        singleBinA.setWaveVectorCoefficients(waveVectorFactoryTarget.getCoefficients());
//        singleBinA.setWaveVectors(waveVectorFactoryTarget.getWaveVectors());
//        singleBinA.setComparedWV(11);
//        System.out.println("singleBinA has set its comparedWV to 11");
//        
//        MeterOverlap singleOverlapinA = new MeterOverlap("singleOverlapinA", Null.DIMENSION, meterAinA, singleBinA, temperature);
//        bMeters[1] = singleOverlapinA;


//REFERENCE
        // Set up REFERENCE system - System B - 0 - Hybrid system
        PotentialMasterList potentialMasterRef = new PotentialMasterList(this, space);
        boxRef = this.makeBox();
        boxRef.setNMolecules(species, numAtoms);

        p2 = new P2HardSphere(space, 1.0, true);
        p2 = new P2XOrder(space, (Potential2HardSpherical) p2);
        p2.setBox(boxRef);
        potentialMasterRef.addPotential(p2, new AtomType[]{
                species.getLeafType(), species.getLeafType()});

        primitive = new PrimitiveCubic(space, 1.0 / density);
        boundaryRef = new BoundaryRectangularPeriodic(space, numAtoms / density);
        nCells = new int[]{numAtoms};
        boxRef.setBoundary(boundaryRef);

        CoordinateDefinitionLeaf coordinateDefinitionRef = new
                CoordinateDefinitionLeaf(boxRef, primitive, space);
        coordinateDefinitionRef.initializeCoordinates(nCells);

        neighborRange = 1.01 / density;
        potentialMasterRef.setRange(neighborRange);
        //find neighbors now.  Don't hook up NeighborListManager since the
        //  neighbors won't change
        potentialMasterRef.getNeighborManager(boxRef).reset();

        IntegratorMC integratorRef = new IntegratorMC(potentialMasterRef,
                random, temperature, boxRef);
        integrators[0] = integratorRef;

        WaveVectorFactory waveVectorFactoryRef = nm.getWaveVectorFactory();
        waveVectorFactoryRef.makeWaveVectors(boxRef);

        compareMove = new MCMoveCompareMultipleWV(potentialMasterRef,
                random);
        integratorRef.getMoveManager().addMCMove(compareMove);
        compareMove.setWaveVectors(waveVectorFactoryRef.getWaveVectors());
        compareMove.setWaveVectorCoefficients(waveVectorFactoryRef.getCoefficients());
        compareMove.setOmegaSquared(nm.getOmegaSquared(),
                waveVectorFactoryRef.getCoefficients());
        compareMove.setEigenVectors(nm.getEigenvectors());
        compareMove.setCoordinateDefinition(coordinateDefinitionRef);
        compareMove.setTemperature(temperature);
        compareMove.setBox(boxRef);
        compareMove.setStepSizeMin(0.001);
        compareMove.setStepSize(0.01);

        meterAinB = new MeterPotentialEnergy(potentialMasterRef, boxRef);

        meterBinB = new MeterCompareMultipleWVBrute(potentialMasterRef,
                coordinateDefinitionRef, boxRef);
        meterBinB.setCoordinateDefinition(coordinateDefinitionRef);
        meterBinB.setEigenVectors(nm.getEigenvectors());
        meterBinB.setOmegaSquared(nm.getOmegaSquared());
        meterBinB.setTemperature(temperature);
        meterBinB.setWaveVectorCoefficients(waveVectorFactoryRef.getCoefficients());
        meterBinB.setWaveVectors(waveVectorFactoryRef.getWaveVectors());
        integratorRef.setMeterPotentialEnergy(meterBinB);

        MeterOverlap meterOverlapInB = new MeterOverlap("MeterOverlapInB", Null.DIMENSION,
                meterBinB, meterAinB, temperature);
        meters[0] = meterOverlapInB;
        potentialMasterRef.getNeighborManager(boxRef).reset();


//JOINT
        //Set up the rest of the joint stuff
        setComparedWV(compWV);
        setChangeableWVs(chbleWV);

        integratorSim = new IntegratorOverlap(new
                IntegratorMC[]{integratorRef, integratorTarget});

        setAccumulator(new AccumulatorVirialOverlapSingleAverage(10, 11, true), 0);
        setAccumulator(new AccumulatorVirialOverlapSingleAverage(10, 11, false), 1);

        setBennettParameter(1.0, 30);

        activityIntegrate = new ActivityIntegrate(integratorSim, 0, true);
        getController().addAction(activityIntegrate);

    }

    public static void main(String args[]) {
        SimOverlapMultipleWaveVectorsParam params = new SimOverlapMultipleWaveVectorsParam();
        String inputFilename = null;
        if (args.length > 0) {
            inputFilename = args[0];
        }
        if (inputFilename != null) {
            ReadParameters readParameters = new
                    ReadParameters(inputFilename, params);
            readParameters.readParameters();
        }

        int numMolecules = params.numAtoms;
        double density = params.density;
        int D = params.D;
        double harmonicFudge = params.harmonicFudge;
        String filename = params.filename;
        if (filename.length() == 0) {
            filename = "1DHR";
        }
        double temperature = params.temperature;
        int[] comparedWV = params.comparedWV;
        int[] changeableWVs = params.changeableWV;

        int numSteps = params.numSteps;
        int runBlockSize = params.runBlockSize;
        int subBlockSize = params.subBlockSize;

        int numEqSteps = params.eqNumSteps;
        int eqBlockSize = params.eqBlockSize;

        int numBenSteps = params.bennettNumSteps;
        int benBlockSize = params.benBlockSize;

        String refFileName = args.length > 0 ? filename + "_ref" : null;


        //instantiate simulations!
        SimOverlapMultipleWaveVectors sim = new SimOverlapMultipleWaveVectors(Space.getInstance(D), numMolecules,
                density, temperature, filename, harmonicFudge, comparedWV, changeableWVs);
        System.out.println(numMolecules + " atoms at density " + density);
        System.out.println("harmonic fudge: " + harmonicFudge);
        System.out.println("temperature: " + temperature);
        System.out.println("Total steps: " + numSteps + " , split into blocks of " + runBlockSize);
        System.out.println(subBlockSize + " steps in subintegrator, per step in  main integrator");
        System.out.println(numEqSteps + " equilibration steps, split into blocks of " + eqBlockSize);
        System.out.println(numBenSteps + " Bennett-only steps, split into blocks of " + benBlockSize);
        System.out.println("output data to " + filename);
        System.out.println("instantiated");

        //Divide out all the steps, so that the subpieces have the proper # of steps
        numSteps /= subBlockSize;
        numEqSteps /= subBlockSize;
        numBenSteps /= subBlockSize;

        //start simulation & equilibrate
        sim.integratorSim.getMoveManager().setEquilibrating(true);
        sim.integratorSim.setNumSubSteps(subBlockSize);

        System.out.println("Init Bennett");
        sim.setAccumulatorBlockSize(benBlockSize);
        sim.initBennettParameter(filename, numBenSteps, benBlockSize);
        if (Double.isNaN(sim.bennettParam) || sim.bennettParam == 0 ||
                Double.isInfinite(sim.bennettParam)) {
            throw new RuntimeException("Simulation failed to find a valid " +
                    "Bennett parameter");
        }

        System.out.println("equilibrate");
        sim.setAccumulatorBlockSize(eqBlockSize);
        sim.equilibrate(refFileName, numEqSteps, eqBlockSize);
        if (Double.isNaN(sim.bennettParam) || sim.bennettParam == 0 ||
                Double.isInfinite(sim.bennettParam)) {
            throw new RuntimeException("Simulation failed to find a valid " +
                    "Bennett parameter");
        }
        System.out.println("equilibration finished.");
//        sim.setBennettParameter(0.573265415766427);
        sim.setAccumulatorBlockSize(runBlockSize);

        sim.integratorSim.getMoveManager().setEquilibrating(false);
        sim.activityIntegrate.setMaxSteps(numSteps);
        sim.getController().actionPerformed();
        System.out.println("final reference optimal step frequency " +
                sim.integratorSim.getIdealRefStepFraction() + " (actual: " +
                sim.integratorSim.getRefStepFraction() + ")");

        //CALCULATION OF HARMONIC ENERGY
        double AHarmonic = CalcHarmonicA.doit(sim.nm, D, temperature, numMolecules);

        double[] ratioAndError = sim.dsvo.getOverlapAverageAndError();
        double ratio = ratioAndError[0];
        double error = ratioAndError[1];
        System.out.println("ratio average: " + ratio + ", error: " + error);
        System.out.println("free energy difference: " + (-Math.log(ratio)) +
                ", error: " + (error / ratio));
        System.out.println("target free energy: " + (AHarmonic - Math.log(ratio)));
        DataGroup allYourBase =
                (DataGroup) sim.accumulators[0].getData(sim.dsvo.minDiffLocation());
        System.out.println("harmonic ratio average: " +
                ((DataDoubleArray) allYourBase.getData(AccumulatorRatioAverageCovariance.RATIO.index)).getData()[1]
                + " error: " +
                ((DataDoubleArray) allYourBase.getData(AccumulatorRatioAverageCovariance.RATIO_ERROR.index)).getData()[1]);

        allYourBase = (DataGroup) sim.accumulators[1].getData(sim.accumulators[0].getNBennetPoints() -
                sim.dsvo.minDiffLocation() - 1);
        System.out.println("target ratio average: " +
                ((DataDoubleArray) allYourBase.getData(AccumulatorRatioAverageCovariance.RATIO.index)).getData()[1]
                + " error: " +
                ((DataDoubleArray) allYourBase.getData(AccumulatorRatioAverageCovariance.RATIO_ERROR.index)).getData()[1]);

        if (D == 1) {
            double AHR = -(numMolecules - 1) * Math.log(numMolecules / density - numMolecules)
                    + SpecialFunctions.lnFactorial(numMolecules);
            System.out.println("Hard-rod free energy: " + AHR);
        }
    }
    
    public void setBennettParameter(double benParamCenter, double span) {
        bennettParam = benParamCenter;
        accumulators[0].setBennetParam(benParamCenter,span);
        accumulators[1].setBennetParam(benParamCenter,span);
    }
    
    public void setBennettParameter(double newBennettParameter) {
        System.out.println("setting ref pref (explicitly) to "+
                newBennettParameter);
        setAccumulator(new AccumulatorVirialOverlapSingleAverage(1,true),0);
        setAccumulator(new AccumulatorVirialOverlapSingleAverage(1,false),1);
        setBennettParameter(newBennettParameter,1);

    }
    
    public void initBennettParameter(String fileName, int initSteps, int initBlockSize) {
        // benParam = -1 indicates we are searching for an appropriate value
        bennettParam = -1.0;
        integratorSim.getMoveManager().setEquilibrating(true);

        if (fileName != null) {
            try {
                FileReader fileReader = new FileReader(fileName);
                BufferedReader bufReader = new BufferedReader(fileReader);
                String benParamString = bufReader.readLine();
                bennettParam = Double.parseDouble(benParamString);
                bufReader.close();
                fileReader.close();
                System.out.println("setting ref pref (from file) to "+bennettParam);
                setAccumulator(new AccumulatorVirialOverlapSingleAverage(1,true),0);
                setAccumulator(new AccumulatorVirialOverlapSingleAverage(1,false),1);
                setBennettParameter(bennettParam,1);
            }
            catch (IOException e) {
                System.out.println("Bennett parameter not from file");
                // file not there, which is ok.
            }
        }

        if (bennettParam == -1) {

//            int oldBlockSize = blockSize;
//            long newBlockSize = initSteps*integratorSim.getNumSubSteps()/1000;
//            //Make sure the new block size is reasonable.
//            if(newBlockSize < 1000){
//                newBlockSize = 1000;
//            }
//            if(newBlockSize > 1000000){
//                newBlockSize = 1000000;
//            }
//            setAccumulatorBlockSize((int)newBlockSize);

            // equilibrate off the lattice to avoid anomolous contributions
            activityIntegrate.setMaxSteps(initSteps);

            getController().actionPerformed();
            getController().reset();

            setAccumulator(new AccumulatorVirialOverlapSingleAverage(initBlockSize,41,true),0);
            setAccumulator(new AccumulatorVirialOverlapSingleAverage(initBlockSize,41,false),1);
            setBennettParameter(1e40,40);
            activityIntegrate.setMaxSteps(initSteps);

            getController().actionPerformed();
            getController().reset();

            int newMinDiffLoc = dsvo.minDiffLocation();
            bennettParam = accumulators[0].getBennetAverage(newMinDiffLoc)
                /accumulators[1].getBennetAverage(newMinDiffLoc);

            double top = accumulators[0].getBennetAverage(newMinDiffLoc);
            System.out.println("top " + top);
            double bottom = accumulators[1].getBennetAverage(newMinDiffLoc);
            System.out.println("bottom " + bottom);

            if (Double.isNaN(bennettParam) || bennettParam == 0 || Double.isInfinite(bennettParam)) {
                throw new RuntimeException("Simulation failed to find a valid ref pref");
            }
            System.out.println("setting ref pref to "+bennettParam);
//            setAccumulatorBlockSize(oldBlockSize);

            setAccumulator(new AccumulatorVirialOverlapSingleAverage(11,true),0);
            setAccumulator(new AccumulatorVirialOverlapSingleAverage(11,false),1);
            setBennettParameter(bennettParam,5);

            // set benParam back to -1 so that later on we know that we've been looking for
            // the appropriate value
            bennettParam = -1;
            getController().reset();
        }
        integratorSim.getMoveManager().setEquilibrating(false);
    }

    public void setAccumulator(AccumulatorVirialOverlapSingleAverage
            newAccumulator, int iBox) {
        accumulators[iBox] = newAccumulator;
        if (accumulatorPumps[iBox] == null) {
            accumulatorPumps[iBox] = new DataPump(meters[iBox], newAccumulator);
            IntegratorListenerAction pumpListener = new IntegratorListenerAction(accumulatorPumps[iBox]);
            pumpListener.setInterval(1);
            integrators[iBox].getEventManager().addListener(pumpListener);
//            integrators[iBox].setActionInterval(accumulatorPumps[iBox],
//                    boxRef.getLeafList().getAtomCount()*2);
        }
        else {
            accumulatorPumps[iBox].setDataSink(newAccumulator);
        }
        if (integratorSim != null && accumulators[0] != null &&
                accumulators[1] != null) {
            dsvo = new DataSourceVirialOverlap(accumulators[0],accumulators[1]);
            integratorSim.setReferenceFracSource(dsvo);
        }

    }

    public void setAccumulatorBlockSize(int newBlockSize) {
        for (int i=0; i<2; i++) {
            accumulators[i].setBlockSize(newBlockSize);
        }
        try {
            // reset the integrator so that it will re-adjust step frequency
            // and ensure it will take enough data for both ref and target
            integratorSim.reset();
        }
        catch (ConfigurationOverlapException e) { /* meaningless */ }
    }

    public void equilibrate(String fileName, int initSteps, int initBlockSize) {
        // run a short simulation to get reasonable MC Move step sizes and
        // (if needed) narrow in on a reference preference
        activityIntegrate.setMaxSteps(initSteps);

        integratorSim.getMoveManager().setEquilibrating(true);

        //This code allows the computer to set the block size for the main
        //simulation and equilibration/finding alpha separately.
//        int oldBlockSize = blockSize;
//        long newBlockSize = initSteps*integratorSim.getNumSubSteps()/1000;
//        //make sure new block size is reasonablel
//        if(newBlockSize < 1000){
//            newBlockSize = 1000;
//        }
//        if (newBlockSize >1000000) {
//            newBlockSize = 1000000;
//        }
//        setAccumulatorBlockSize((int)newBlockSize);

//        setAccumulatorBlockSize((int)eqBlockSize);

        for (int i=0; i<2; i++) {
            if (integrators[i] instanceof IntegratorMC) integrators[i].getMoveManager().setEquilibrating(true);
        }
        getController().actionPerformed();
        getController().reset();
        for (int i=0; i<2; i++) {
            if (integrators[i] instanceof IntegratorMC) integrators[i].getMoveManager().setEquilibrating(false);
        }

        if (bennettParam == -1) {
            int newMinDiffLoc = dsvo.minDiffLocation();
            bennettParam = accumulators[0].getBennetAverage(newMinDiffLoc)
                /accumulators[1].getBennetAverage(newMinDiffLoc);
            System.out.println("setting ref pref to "+bennettParam+" ("+newMinDiffLoc+")");
            setAccumulator(new AccumulatorVirialOverlapSingleAverage(initBlockSize,1,true),0);
            setAccumulator(new AccumulatorVirialOverlapSingleAverage(initBlockSize,1,false),1);
            setBennettParameter(bennettParam,1);
            if (fileName != null) {
                try {
                    FileWriter fileWriter = new FileWriter(fileName);
                    BufferedWriter bufWriter = new BufferedWriter(fileWriter);
                    bufWriter.write(String.valueOf(bennettParam)+"\n");
                    bufWriter.close();
                    fileWriter.close();
                }
                catch (IOException e) {
                    throw new RuntimeException("couldn't write to Bennet parameter file");
                }
            }
        }
        else {
            dsvo.reset();
        }

        integratorSim.getMoveManager().setEquilibrating(false);
//        setAccumulatorBlockSize(oldBlockSize);
    }
    
    public void setComparedWV(int[] cwvs){
        meterBinB.setComparedWV(cwvs);
        meterBinA.setComparedWV(cwvs);
        compareMove.setComparedWV(cwvs);
        changeMove.addChangeableWV(cwvs);
        
        System.out.println("Compared WV: ");
        for (int i = 0; i < cwvs.length; i++){
            System.out.println(cwvs[i]);
        }
    }
    
    public void setChangeableWVs(int[] cwvs){
        changeMove.addChangeableWV(cwvs);
        compareMove.setChangeableWVs(cwvs);
        
        System.out.println("Hard Rod WV: ");
        for (int i = 0; i < cwvs.length; i++){
            System.out.println(cwvs[i]);
        }
    }
    
    
    public static class SimOverlapMultipleWaveVectorsParam extends ParameterBase {
        public int numAtoms = 32;
        public double density = 0.5;
        public int D = 1;
        public double harmonicFudge = 1.0;
        public String filename = "HR1D_";
        public double temperature = 1.0;
        public int[] comparedWV = {2};
        public int[] changeableWV = {1};
        
        public int numSteps = 4000000;
        public int runBlockSize = 100000;
        public int subBlockSize = 1000;    //# of steps in subintegrator per integrator step
        
        public int eqNumSteps = 40000;  
        public int eqBlockSize = 1000;
        
        public int bennettNumSteps = 40000;
        public int benBlockSize = 1000;
    }
}
