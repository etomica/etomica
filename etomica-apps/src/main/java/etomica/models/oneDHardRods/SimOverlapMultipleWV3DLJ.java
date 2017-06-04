/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.models.oneDHardRods;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;

import etomica.action.activity.ActivityIntegrate;
import etomica.atom.IAtomType;
import etomica.box.Box;
import etomica.data.DataPump;
import etomica.data.IEtomicaDataSource;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.data.types.DataDoubleArray;
import etomica.data.types.DataGroup;
import etomica.exception.ConfigurationOverlapException;
import etomica.integrator.IntegratorMC;
import etomica.lattice.crystal.Basis;
import etomica.lattice.crystal.BasisCubicFcc;
import etomica.lattice.crystal.BasisMonatomic;
import etomica.lattice.crystal.Primitive;
import etomica.lattice.crystal.PrimitiveCubic;
import etomica.listener.IntegratorListenerAction;
import etomica.math.SpecialFunctions;
import etomica.normalmode.CoordinateDefinitionLeaf;
import etomica.normalmode.MeterNormalMode;
import etomica.normalmode.NormalModesFromFile;
import etomica.normalmode.WaveVectorFactory;
import etomica.normalmode.WriteS;
import etomica.overlap.IntegratorOverlap;
import etomica.potential.P2LennardJones;
import etomica.potential.P2SoftSphericalTruncatedShifted;
import etomica.potential.Potential2SoftSpherical;
import etomica.potential.PotentialMasterMonatomic;
import etomica.simulation.Simulation;
import etomica.space.Boundary;
import etomica.space.BoundaryDeformableLattice;
import etomica.space.Space;
import etomica.species.SpeciesSpheresMono;
import etomica.units.Null;
import etomica.util.ParameterBase;
import etomica.util.ReadParameters;
import etomica.virial.overlap.AccumulatorVirialOverlapSingleAverage;
import etomica.virial.overlap.DataSourceVirialOverlap;

public class SimOverlapMultipleWV3DLJ extends Simulation {

    private static final String APP_NAME = "SimOverlapMultipleWV3DLJ";
    Primitive primitiveTarget, primitiveRef;
    int[] nCells;
    NormalModesFromFile nm;
    double bennettParam;       //adjustable parameter - Bennett's parameter
    public IntegratorOverlap integratorSim; //integrator for the whole simulation
    public DataSourceVirialOverlap dsvo;
    public BasisMonatomic basis;
    ActivityIntegrate activityIntegrate;
    
    IntegratorMC[] integrators;
    public AccumulatorVirialOverlapSingleAverage[] accumulators;
    public DataPump[] accumulatorPumps;
    public IEtomicaDataSource[] meters;
    public Box boxTarget, boxRef;
    public Boundary boundaryTarget, boundaryRef;
    MCMoveChangeMultipleWV changeMove;
    MCMoveCompareMultipleWV compareMove;
    MeterPotentialEnergy meterAinB, meterAinA;
    MeterCompareMultipleWVBrute meterBinA, meterBinB;
    
    MeterNormalMode mnm;
    WriteS sWriter;
    
    public SimOverlapMultipleWV3DLJ(Space _space, int numAtoms, double 
            density, double temperature, String filename, double harmonicFudge,
            int[] compWV, int[] harmWV){
        super(_space);
        
        System.out.println("THIS CODE IS NOT FINISHED!");
        System.out.println("need to fix this setHarmonicWV");
        
        
//        long seed = 2;
//        System.out.println("Seed explicitly set to " + seed);
//        IRandom rand = new RandomNumberGenerator(seed);
//        this.setRandom(rand);
        
        //Set up some of the joint stuff
        SpeciesSpheresMono species = new SpeciesSpheresMono(this, space);
        addSpecies(species);
        
        integrators = new IntegratorMC[2];
        accumulatorPumps = new DataPump[2];
        meters = new IEtomicaDataSource[2];
        accumulators = new AccumulatorVirialOverlapSingleAverage[2];
        
        basis = new BasisMonatomic(space);
        
//TARGET    
        // Set up target system - A, 1, hard rod
        PotentialMasterMonatomic potentialMasterTarget = new 
                PotentialMasterMonatomic(this);
        boxTarget = new Box(space);
        addBox(boxTarget);
        boxTarget.setNMolecules(species, numAtoms);
        
        primitiveTarget = new PrimitiveCubic(space, 1.0);
        double v = primitiveTarget.unitCell().getVolume();
        primitiveTarget.scaleSize(Math.pow(v*density/4, -1.0/3.0));
        int numberOfCells = (int)Math.round(Math.pow(numAtoms/4, 1.0/3.0));
        nCells = new int[]{numberOfCells, numberOfCells, numberOfCells};
        boundaryTarget = new BoundaryDeformableLattice(primitiveTarget, nCells);
        boxTarget.setBoundary(boundaryTarget);
        Basis basisTarget = new BasisCubicFcc();
        
        CoordinateDefinitionLeaf coordinateDefinitionTarget = new 
                CoordinateDefinitionLeaf(boxTarget, primitiveTarget, basisTarget, space);
        coordinateDefinitionTarget.initializeCoordinates(nCells);
        
        Potential2SoftSpherical p2 = new P2LennardJones(space, 1.0, 1.0);
        double truncationRadius = boundaryTarget.getBoxSize().getX(0) * 0.495;
        P2SoftSphericalTruncatedShifted pTruncated = new 
                P2SoftSphericalTruncatedShifted(space, p2, truncationRadius);
        potentialMasterTarget.addPotential(pTruncated, new IAtomType[]
                {species.getLeafType(), species.getLeafType()});
        
        IntegratorMC integratorTarget = new IntegratorMC(potentialMasterTarget,
                random, temperature);
        integrators[1] = integratorTarget;
        integratorTarget.setBox(boxTarget);
        
        nm = new NormalModesFromFile(filename, space.D());
        nm.setHarmonicFudge(harmonicFudge);
        nm.setTemperature(temperature);
        nm.getOmegaSquared();
        
        WaveVectorFactory waveVectorFactoryTarget = nm.getWaveVectorFactory();
        waveVectorFactoryTarget.makeWaveVectors(boxTarget);
        int wvflength = waveVectorFactoryTarget.getWaveVectors().length;
        System.out.println("We have " + wvflength +" wave vectors.");
        System.out.println("Wave Vector Coefficients:");
        for(int i = 0; i < wvflength; i++){
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
        
        meterAinA = new MeterPotentialEnergy(potentialMasterTarget);
        meterAinA.setBox(boxTarget);
        
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
        PotentialMasterMonatomic potentialMasterRef = new 
                PotentialMasterMonatomic(this);
        boxRef = new Box(space);
        addBox(boxRef);
        boxRef.setNMolecules(species, numAtoms);
        
        primitiveRef = new PrimitiveCubic(space, 1.0);
        v = primitiveRef.unitCell().getVolume();
        primitiveRef.scaleSize(Math.pow(v*density/4, -1.0/3.0));
        numberOfCells = (int)Math.round(Math.pow(numAtoms/4, 1.0/3.0));
        nCells = new int[]{numberOfCells, numberOfCells, numberOfCells};
        boundaryRef  = new BoundaryDeformableLattice(primitiveRef, nCells);
        boxRef.setBoundary(boundaryRef);
        Basis basisRef = new BasisCubicFcc();
        
        CoordinateDefinitionLeaf coordinateDefinitionRef = new 
                CoordinateDefinitionLeaf(boxRef, primitiveRef, space);
        coordinateDefinitionRef.initializeCoordinates(nCells);
        
        p2 = new P2LennardJones(space, 1.0, 1.0);
        truncationRadius = boundaryTarget.getBoxSize().getX(0) * 0.5;
        pTruncated = new P2SoftSphericalTruncatedShifted(space, p2, truncationRadius);
        potentialMasterRef.addPotential(pTruncated, new IAtomType[]
                {species.getLeafType(), species.getLeafType()});
        
        IntegratorMC integratorRef = new IntegratorMC(potentialMasterRef, 
                random, temperature);
        integratorRef.setBox(boxRef);
        integrators[0] = integratorRef;
        
        nm = new NormalModesFromFile(filename, space.D());
        nm.setHarmonicFudge(harmonicFudge);
        nm.setTemperature(temperature);
        
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
        
        meterAinB = new MeterPotentialEnergy(potentialMasterRef);
        meterAinB.setBox(boxRef);
       
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
        
        integratorRef.setBox(boxRef);
        
        //Stuff to take care of recording spring constant!!
        mnm = new MeterNormalMode();
        mnm.setCoordinateDefinition(coordinateDefinitionRef);
        mnm.setWaveVectorFactory(waveVectorFactoryRef);
        mnm.setBox(boxRef);
        
        IntegratorListenerAction mnmListener = new IntegratorListenerAction(mnm);
        integratorRef.getEventManager().addListener(mnmListener);
        mnmListener.setInterval(1000);
        
        sWriter = new WriteS(space);
        sWriter.setFilename(filename + "_output_" + compWV[0]);
        sWriter.setMeter(mnm);
        sWriter.setWaveVectorFactory(mnm.getWaveVectorFactory());
        sWriter.setOverwrite(true);
        
//        meterBinB.getSingle().setCoordinateDefinition(coordinateDefinitionRef);
//        meterBinB.getSingle().setEigenVectors(nm.getEigenvectors(boxRef));
//        meterBinB.getSingle().setOmegaSquared(nm.getOmegaSquared(boxRef));
//        meterBinB.getSingle().setTemperature(temperature);
//        meterBinB.getSingle().setWaveVectorCoefficients(waveVectorFactoryRef.getCoefficients());
//        meterBinB.getSingle().setWaveVectors(waveVectorFactoryRef.getWaveVectors());
//        meterBinB.setA(false);
        
//        singleBinB = new MeterCompareSingleModeBrute(potentialMasterRef, 
//                coordinateDefinitionRef, boxRef);
//        singleBinB.setCoordinateDefinition(coordinateDefinitionRef);
//        singleBinB.setEigenVectors(nm.getEigenvectors(boxRef));
//        singleBinB.setOmegaSquared(nm.getOmegaSquared(boxRef));
//        singleBinB.setTemperature(temperature);
//        singleBinB.setWaveVectorCoefficients(waveVectorFactoryRef.getCoefficients());
//        singleBinB.setWaveVectors(waveVectorFactoryRef.getWaveVectors());
//        
//        
//        MeterOverlap singleOverlapInB = new MeterOverlap("SingleOverlapInB", Null.DIMENSION, singleBinB, meterAinB, temperature);
//        bMeters[0] = singleOverlapInB;
        
        
        
//JOINT
        //Set up the rest of the joint stuff
        setComparedWV(compWV);
        setHarmonicWV(harmWV);
        
        integratorSim = new IntegratorOverlap(new 
                IntegratorMC[]{integratorRef, integratorTarget});
        
        setAccumulator(new AccumulatorVirialOverlapSingleAverage(10, 11, true), 0);
        setAccumulator(new AccumulatorVirialOverlapSingleAverage(10, 11, false), 1);
        
        setBennettParameter(1.0, 30);
        
        activityIntegrate = new ActivityIntegrate(integratorSim, 0, true);
        getController().addAction(activityIntegrate);
        
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
    
    public static void main(String args[]){
        SimOverlapMultipleWaveVectorsParam params = new SimOverlapMultipleWaveVectorsParam();
        String inputFilename = null;
        if(args.length > 0){
            inputFilename = args[0];
        }
        if(inputFilename != null){
            ReadParameters readParameters = new 
                ReadParameters(inputFilename, params);
            readParameters.readParameters();
        }
        
        int numMolecules = params.numAtoms;
        double density = params.density;
        int D = params.D;
        double harmonicFudge = params.harmonicFudge;
        String filename = params.filename;
        if(filename.length() == 0){
            filename = "1DHR";
        }
        double temperature = params.temperature;
        int[] comparedWV = params.comparedWV;
        int[] harmonicWV = params.harmonicWV;
        
        int numSteps = params.numSteps;
        int runBlockSize = params.runBlockSize;
        int subBlockSize = params.subBlockSize;
        
        int numEqSteps = params.eqNumSteps;
        int eqBlockSize = params.eqBlockSize;
    
        int numBenSteps = params.bennettNumSteps;
        int benBlockSize = params.benBlockSize;
        
        String refFileName = args.length > 0 ? filename+"_ref" : null;
        
        //instantiate simulations!
        SimOverlapMultipleWV3DLJ sim = new SimOverlapMultipleWV3DLJ(Space.getInstance(D), numMolecules,
                density, temperature, filename, harmonicFudge, comparedWV, harmonicWV);
        System.out.println("Running " + sim.APP_NAME);
        System.out.println(numMolecules+" atoms at density "+density);
        System.out.println("harmonic fudge: "+harmonicFudge);
        System.out.println("temperature: " + temperature);
        System.out.println("compared wave vectors: ");
        for(int i = 0; i < comparedWV.length; i++){
            System.out.println(comparedWV[i]);
        }
        System.out.println("harmonic wave vectors: ");
        for(int i = 0; i < harmonicWV.length; i++ ){
            System.out.println(harmonicWV[i]);
        }
        System.out.println("Total steps: "+numSteps+" , split into blocks of "+runBlockSize);
        System.out.println(subBlockSize+" steps in subintegrator, per step in  main integrator");
        System.out.println(numEqSteps+" equilibration steps, split into blocks of "+ eqBlockSize);
        System.out.println(numBenSteps+" Bennett-only steps, split into blocks of "+benBlockSize);
        System.out.println("output data to "+filename);
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
        if(Double.isNaN(sim.bennettParam) || sim.bennettParam == 0 || 
                Double.isInfinite(sim.bennettParam)){
            throw new RuntimeException("Simulation failed to find a valid " +
                    "Bennett parameter");
        }
        
        System.out.println("equilibrate");
        sim.setAccumulatorBlockSize(eqBlockSize);
        sim.equilibrate(refFileName, numEqSteps, eqBlockSize);
        if(Double.isNaN(sim.bennettParam) || sim.bennettParam == 0 || 
                Double.isInfinite(sim.bennettParam)){
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
        
        double[][] omega2 = sim.nm.getOmegaSquared(); 
        //Above known from the analytical results. - otherwise it would be from 
        //the S matrix.
        double[] coeffs = sim.nm.getWaveVectorFactory().getCoefficients();
        
        //CALCULATION OF HARMONIC ENERGY
        double AHarmonic = 0;
        for(int i=0; i<omega2.length; i++) {
            for(int j=0; j<omega2[0].length; j++) {
                if (!Double.isInfinite(omega2[i][j])) {
                    AHarmonic += coeffs[i] * Math.log(omega2[i][j]*coeffs[i] /
                            (temperature*Math.PI));
                }
            }
        }
        int totalCells = 1;
        for (int i=0; i<D; i++) {
            totalCells *= sim.nCells[i];
        }
        int basisSize = sim.basis.getScaledCoordinates().length;
        double fac = 1;
        if (totalCells % 2 == 0) {
            fac = Math.pow(2,D);
        }
        AHarmonic -= Math.log(Math.pow(2.0, basisSize *D * (totalCells - fac) / 
                2.0) / Math.pow(totalCells, 0.5 * D));
        System.out.println("Harmonic-reference free energy: " + AHarmonic * 
                temperature);
        double[] ratioAndError = sim.dsvo.getOverlapAverageAndError();
        double ratio = ratioAndError[0];
        double error = ratioAndError[1];
        System.out.println("ratio average: "+ratio+", error: "+error);
        System.out.println("free energy difference: " + (-Math.log(ratio)) + 
                ", error: "+(error/ratio));
        System.out.println("target free energy: " + (AHarmonic-Math.log(ratio)));
        DataGroup allYourBase = 
            (DataGroup)sim.accumulators[0].getData(sim.dsvo.minDiffLocation());
        System.out.println("harmonic ratio average: " + 
                ((DataDoubleArray)allYourBase.getData(sim.accumulators[0].RATIO.index)).getData()[1]
                 + " error: " + 
                ((DataDoubleArray)allYourBase.getData(sim.accumulators[0].RATIO_ERROR.index)).getData()[1]);
        
        allYourBase = (DataGroup)sim.accumulators[1].getData(sim.accumulators[1].getNBennetPoints() -
                sim.dsvo.minDiffLocation()-1);
        System.out.println("target ratio average: " + 
                ((DataDoubleArray)allYourBase.getData(sim.accumulators[1].RATIO.index)).getData()[1]
                 + " error: " + 
                ((DataDoubleArray)allYourBase.getData(sim.accumulators[1].RATIO_ERROR.index)).getData()[1]);
    
        if(D==1) {
            double AHR = -(numMolecules-1)*Math.log(numMolecules/density-numMolecules)
                + SpecialFunctions.lnFactorial(numMolecules) ;
            System.out.println("Hard-rod free energy: "+AHR);
        }
    }
    
    public void setComparedWV(int[] cwvs){
        meterBinB.setComparedWV(cwvs);
        meterBinA.setComparedWV(cwvs);
        compareMove.setComparedWV(cwvs);
        
//        meterBinB.getSingle().setComparedWV(cwvs[0]);
//        meterBinA.getSingle().setComparedWV(cwvs[0]);
        
        
    }
    
    public void setHarmonicWV(int[] hwv){
        System.out.println("THIS CODE IS NOT FINISHED!");
        System.out.println("need to fix this setHarmonicWV");
        
        
//        compareMove.setHarmonicWV(hwv);
//        changeMove.setHarmonicWV(hwv);
    }
    
    public static class SimOverlapMultipleWaveVectorsParam extends ParameterBase {
        public int numAtoms = 10;
        public double density = 1.3;
        public int D = 3;
        public double harmonicFudge = 1.0;
        public String filename = "normal_modes_LJ_3D_32";
        public double temperature = 1.0;
        public int[] comparedWV = {1, 2};
        public int[] harmonicWV = {3, 4, 5, 6, 7};
        
        public int numSteps = 400000;
        public int runBlockSize = 1000;
        public int subBlockSize = 1000;    //# of steps in subintegrator per integrator step
        
        public int eqNumSteps = 40000;  
        public int eqBlockSize = 1000;
        
        public int bennettNumSteps = 40000;
        public int benBlockSize = 1000;
    }
}
