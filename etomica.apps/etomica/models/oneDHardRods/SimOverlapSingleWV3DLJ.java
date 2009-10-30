package etomica.models.oneDHardRods;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;

import etomica.action.activity.ActivityIntegrate;
import etomica.api.IAtomType;
import etomica.api.IBox;
import etomica.box.Box;
import etomica.data.AccumulatorRatioAverage;
import etomica.data.DataPump;
import etomica.data.DataSourceScalar;
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
import etomica.normalmode.NormalModes;
import etomica.normalmode.NormalModesFromFile;
import etomica.normalmode.WaveVectorFactory;
import etomica.normalmode.WriteS;
import etomica.potential.P2LennardJones;
import etomica.potential.P2SoftSphericalTruncatedShifted;
import etomica.potential.Potential2SoftSpherical;
import etomica.potential.PotentialMasterMonatomic;
import etomica.simulation.Simulation;
import etomica.space.Boundary;
import etomica.space.BoundaryDeformableLattice;
import etomica.space.Space;
import etomica.species.SpeciesSpheresMono;
import etomica.units.Energy;
import etomica.units.Null;
import etomica.util.ParameterBase;
import etomica.util.ReadParameters;
import etomica.virial.overlap.AccumulatorVirialOverlapSingleAverage;
import etomica.virial.overlap.DataSourceVirialOverlap;
import etomica.virial.overlap.IntegratorOverlap;

/**
 * MC simulation of Lennard-Jones system.
 * @author cribbin
 */
public class SimOverlapSingleWV3DLJ extends Simulation {
    private static final long serialVersionUID = 1L;
    private static final String APP_NAME = "SimSingleWaveVector";
    Primitive primitiveTarget, primitiveRef;
    int[] nCells;
    NormalModes nm;
    double bennettParam;       //adjustable parameter - Bennett's parameter
    public IntegratorOverlap integratorSim; //integrator for the whole simulation
    public DataSourceVirialOverlap dsvo;
    public BasisMonatomic basis;
    ActivityIntegrate activityIntegrate;
    
    IntegratorMC[] integrators;
    
    public AccumulatorVirialOverlapSingleAverage[] accumulators;
    public DataPump[] accumulatorPumps;
    public IEtomicaDataSource[] meters;
    public IBox boxTarget, boxRef;
    public Boundary boundaryTarget, boundaryRef;
    MCMoveChangeSingleWVLoop changeMove;
    MCMoveCompareSingleWVLoop compareMove;
    MeterPotentialEnergy meterAinB, meterAinA;
    MeterCompareSingleWVBrute meterBinA, meterBinB;
    MeterPotentialEnergyDifference meterWrapAinA, meterWrapAinB, meterWrapBinA, meterWrapBinB;
    
    MeterNormalMode mnm;
    WriteS sWriter;
    
    public SimOverlapSingleWV3DLJ(Space _space, int numAtoms, double density, double 
            temperature, String filename, double harmonicFudge, int awv){
        super(_space);
        
//        long seed = 5;
//        System.out.println("Seed explicitly set to " + seed);
//        IRandom rand = new RandomNumberGenerator(seed);
//        this.setRandom(rand);
        
        //Set up some of the joint stuff
        SpeciesSpheresMono species = new SpeciesSpheresMono(this, space);
        getSpeciesManager().addSpecies(species);
        
        integrators = new IntegratorMC[2];
        accumulatorPumps = new DataPump[2];
        meters = new IEtomicaDataSource[2];
        accumulators = new AccumulatorVirialOverlapSingleAverage[2];
        
        basis = new BasisMonatomic(space);
        
//TARGET
        //Set up target system   - A - 1
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
        
        CoordinateDefinitionLeaf coordinateDefinitionTarget = new CoordinateDefinitionLeaf(boxTarget, primitiveTarget, basisTarget, space);
        coordinateDefinitionTarget.initializeCoordinates(nCells);
        
//        for(int k = 0; k < 32; k++){
//            System.out.println(k + " " +((IAtomPositioned)coordinateDefinitionTarget.getBox().getLeafList().getAtom(k)).getPosition());
//        }
        
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
        
        changeMove = new MCMoveChangeSingleWVLoop(potentialMasterTarget, random);
        integratorTarget.getMoveManager().addMCMove(changeMove);
        changeMove.setWaveVectors(waveVectorFactoryTarget.getWaveVectors());
        changeMove.setWaveVectorCoefficients(
                waveVectorFactoryTarget.getCoefficients());
        changeMove.setEigenVectors(nm.getEigenvectors());
        changeMove.setCoordinateDefinition(coordinateDefinitionTarget);
        changeMove.setBox((IBox)boxTarget);
        changeMove.setStepSizeMin(0.001);
        changeMove.setStepSize(0.01);
        changeMove.setOmegaSquared(nm.getOmegaSquared());
        
        meterAinA = new MeterPotentialEnergy(potentialMasterTarget);
        meterAinA.setBox(boxTarget);
        double latticeEnergy = meterAinA.getDataAsScalar();
        meterWrapAinA = new MeterPotentialEnergyDifference(meterAinA, latticeEnergy);
        
        meterBinA = new MeterCompareSingleWVBrute("meterBinA", potentialMasterTarget, 
                coordinateDefinitionTarget, boxTarget);
        meterBinA.setEigenVectors(nm.getEigenvectors());
        meterBinA.setOmegaSquared(nm.getOmegaSquared());
        meterBinA.setTemperature(temperature);
        meterBinA.setWaveVectorCoefficients(waveVectorFactoryTarget.getCoefficients());
        meterBinA.setWaveVectors(waveVectorFactoryTarget.getWaveVectors());
        meterWrapBinA = new MeterPotentialEnergyDifference(meterBinA, latticeEnergy);
        
        MeterOverlap meterOverlapInA = new MeterOverlap("MeterOverlapInA", Null.DIMENSION, 
                meterWrapAinA, meterWrapBinA, temperature);
        meters[1] = meterOverlapInA;
        
        
//REFERENCE        
        //Set up REFERENCE system - System B - 0 - hybrid system
        PotentialMasterMonatomic potentialMasterRef = new 
            PotentialMasterMonatomic(this);
        boxRef = new Box(space);
        addBox(boxRef);
        boxRef.setNMolecules(species, numAtoms);
//        accumulators[1] = new AccumulatorVirialOverlapSingleAverage(10, 11, true);

        primitiveRef = new PrimitiveCubic(space, 1.0);
        v = primitiveRef.unitCell().getVolume();
        primitiveRef.scaleSize(Math.pow(v*density/4, -1.0/3.0));
        numberOfCells = (int)Math.round(Math.pow(numAtoms/4, 1.0/3.0));
        nCells = new int[]{numberOfCells, numberOfCells, numberOfCells};
        boundaryRef  = new BoundaryDeformableLattice(primitiveRef, nCells);
        boxRef.setBoundary(boundaryRef);
        Basis basisRef = new BasisCubicFcc();
        
        CoordinateDefinitionLeaf coordinateDefinitionRef = new 
                CoordinateDefinitionLeaf(boxRef, primitiveRef, basisRef, space);
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
        
        compareMove = new MCMoveCompareSingleWVLoop(potentialMasterRef, 
                random);
        integratorRef.getMoveManager().addMCMove(compareMove);
        compareMove.setWaveVectors(waveVectorFactoryRef.getWaveVectors());
        compareMove.setWaveVectorCoefficients(waveVectorFactoryRef.getCoefficients());
        compareMove.setOmegaSquared(nm.getOmegaSquared(), 
                waveVectorFactoryRef.getCoefficients());
        compareMove.setEigenVectors(nm.getEigenvectors());
        compareMove.setCoordinateDefinition(coordinateDefinitionRef);
        compareMove.setTemperature(temperature);
        compareMove.setBox((IBox)boxRef);
        compareMove.setStepSizeMin(0.001);
        compareMove.setStepSize(0.01);
        compareMove.setOmegaSquared(nm.getOmegaSquared(), 
                nm.getWaveVectorFactory().getCoefficients());
        
        meterAinB = new MeterPotentialEnergy(potentialMasterRef);
        meterAinB.setBox(boxRef);
        meterWrapAinB = new MeterPotentialEnergyDifference(meterAinB, latticeEnergy);
       
        meterBinB = new MeterCompareSingleWVBrute(potentialMasterRef,
                coordinateDefinitionRef, boxRef);
        meterBinB.setCoordinateDefinition(coordinateDefinitionRef);
        meterBinB.setEigenVectors(nm.getEigenvectors());
        meterBinB.setOmegaSquared(nm.getOmegaSquared());
        meterBinB.setTemperature(temperature);
        meterBinB.setWaveVectorCoefficients(waveVectorFactoryRef.getCoefficients());
        meterBinB.setWaveVectors(waveVectorFactoryRef.getWaveVectors());
        integratorRef.setMeterPotentialEnergy(meterBinB);
        meterWrapBinB = new MeterPotentialEnergyDifference(meterBinB, latticeEnergy);
        
        MeterOverlap meterOverlapInB = new MeterOverlap("MeterOverlapInB", 
                Null.DIMENSION, meterWrapBinB, meterWrapAinB, temperature);
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
        sWriter.setFilename(filename + "_output_" + awv);
        sWriter.setMeter(mnm);
        sWriter.setWaveVectorFactory(mnm.getWaveVectorFactory());
        sWriter.setOverwrite(true);
        
//JOINT
        //Set up the rest of the joint stuff
        setComparedWV(awv);
        
        integratorSim = new IntegratorOverlap(new 
                IntegratorMC[]{integratorRef, integratorTarget});
        
        setAccumulator(new AccumulatorVirialOverlapSingleAverage(10, 11, true), 0);
        setAccumulator(new AccumulatorVirialOverlapSingleAverage(10, 11, false), 1);
        
        setBennettParameter(1.0, 30);
        
        activityIntegrate = new ActivityIntegrate(integratorSim, 0, true);
        getController().addAction(activityIntegrate);
        
    }
    
    /*
     * Sets you muliple bennett parameters.
     */
    public void setBennettParameter(double benParamCenter, double span) {
        bennettParam = benParamCenter;
        accumulators[0].setBennetParam(benParamCenter,span);
        accumulators[1].setBennetParam(benParamCenter,span);
    }
    
    /*
     * Sets you a single bennett parameter
     */
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
            setAccumulatorBlockSize(initBlockSize);
            
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
            if (Double.isNaN(bennettParam) || bennettParam == 0 || Double.isInfinite(bennettParam)) {
                if(Double.isNaN(bennettParam)){ System.out.println("isNan"); }
                if(bennettParam == 0){ System.out.println("zero"); };
                if(Double.isInfinite(bennettParam)){ System.out.println("isinfinite"); }
                
                System.out.println("numerator " + accumulators[0].getBennetAverage(newMinDiffLoc));
                System.out.println("denominator " + accumulators[1].getBennetAverage(newMinDiffLoc));
                
                throw new RuntimeException("Simulation failed to find a valid ref pref");
            }
            System.out.println("setting ref pref to "+bennettParam);
//            setAccumulatorBlockSize(oldBlockSize);
            
            setAccumulator(new AccumulatorVirialOverlapSingleAverage(initBlockSize,11,true),0);
            setAccumulator(new AccumulatorVirialOverlapSingleAverage(initBlockSize,11,false),1);
            setBennettParameter(bennettParam,5);

            // set benParam back to -1 so that later on we know that we've been looking for
            // the appropriate value
            bennettParam = -1;
            getController().reset();
        }
        integratorSim.getMoveManager().setEquilibrating(false);
//        setAccumulatorBlockSize(runBlockSize);
//        activityIntegrate.setMaxSteps(numSteps);
    }
    public void setAccumulator(AccumulatorVirialOverlapSingleAverage 
            newAccumulator, int iBox) {
        accumulators[iBox] = newAccumulator;
//        accumulators[iBox].setBlockSize(blockSize);
//        System.out.println("setAccumlator set to " + blockSize + " blocksize");
        if (accumulatorPumps[iBox] == null) {
            accumulatorPumps[iBox] = new DataPump(meters[iBox], newAccumulator);
            IntegratorListenerAction pumpListener = new IntegratorListenerAction(accumulatorPumps[iBox]);
            integrators[iBox].getEventManager().addListener(pumpListener);
//            integrators[iBox].setActionInterval(accumulatorPumps[iBox], 
//                    boxRef.getLeafList().getAtomCount()*2);
            pumpListener.setInterval(1);
        }
        else {
            accumulatorPumps[iBox].setDataSink(newAccumulator);
        }
        if (integratorSim != null && accumulators[0] != null && 
                accumulators[1] != null) {
            dsvo = new DataSourceVirialOverlap(accumulators[0],accumulators[1]);
            integratorSim.setDSVO(dsvo);
        }
    }
    
    public void setAccumulatorBlockSize(int newBlockSize) {
//        blockSize = newBlockSize;
        for (int i=0; i<2; i++) {
            accumulators[i].setBlockSize(newBlockSize);
//            System.out.println("setAccumlatorBlockSize [] set to " + newBlockSize + " blocksize");
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
        setAccumulatorBlockSize(initBlockSize);
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
//        setAccumulatorBl ockSize((int)newBlockSize);
        
        for (int i=0; i<2; i++) {
            if (integrators[i] instanceof IntegratorMC) ((IntegratorMC)integrators[i]).getMoveManager().setEquilibrating(true);
        }
        getController().actionPerformed();
        getController().reset();
        for (int i=0; i<2; i++) {
            if (integrators[i] instanceof IntegratorMC) ((IntegratorMC)integrators[i]).getMoveManager().setEquilibrating(false);
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
    }
    
    public static void main(String args[]){
        SimOverlapSingleWaveVector3DParam params = new SimOverlapSingleWaveVector3DParam();
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
        int comparedWV = params.comparedWV;
        
        String refFileName = args.length > 0 ? filename+"_ref" : null;
        
        //instantiate simulations!
        SimOverlapSingleWV3DLJ sim = new SimOverlapSingleWV3DLJ  (Space.getInstance(D), numMolecules,
                density, temperature, filename, harmonicFudge, comparedWV);
        int numSteps = params.numSteps;
        int runBlockSize = params.runBlockSize;
        int subBlockSize = params.subBlockSize;
        
        int numEqSteps = params.eqNumSteps;
        int eqBlockSize = params.eqBlockSize;
    
        int numBenSteps = params.bennettNumSteps;
        int benBlockSize = params.benBlockSize;
        
        System.out.println("Running Nancy's single " +D+"D Lennard Jones simulation");
        System.out.println(numMolecules+" atoms at density "+density);
        System.out.println("harmonic fudge: "+harmonicFudge);
        System.out.println("temperature: " + temperature);
        System.out.println("compared wave vector: " + comparedWV);
        System.out.println("Total steps: "+params.numSteps+" , split into blocks of "+ runBlockSize);
        System.out.println(params.subBlockSize+" steps in subintegrator, per step in  main integrator");
        System.out.println(params.eqNumSteps+" equilibration steps, split into blocks of "+ params.eqBlockSize);
        System.out.println(params.bennettNumSteps +" Bennett-only steps, split into blocks of "+params.benBlockSize);
        System.out.println("output data to "+filename);

        System.out.println("instantiated");
        
        //Divide out all the steps, so that the subpieces have the proper # of steps
        numSteps /= params.subBlockSize;
        numEqSteps /= params.subBlockSize;
        numBenSteps /= params.subBlockSize;
        
        //start simulation & equilibrate
        sim.integratorSim.getMoveManager().setEquilibrating(true);
        sim.integratorSim.setNumSubSteps(subBlockSize);
        
        sim.setAccumulatorBlockSize(benBlockSize);
        sim.initBennettParameter(filename, numBenSteps, benBlockSize);
        if(Double.isNaN(sim.bennettParam) || sim.bennettParam == 0 || 
                Double.isInfinite(sim.bennettParam)){
            throw new RuntimeException("Simulation failed to find a valid " +
                    "Bennett parameter");
        }
        sim.setAccumulatorBlockSize(eqBlockSize);
        sim.equilibrate(refFileName, numEqSteps, eqBlockSize);
        if(Double.isNaN(sim.bennettParam) || sim.bennettParam == 0 || 
                Double.isInfinite(sim.bennettParam)){
            throw new RuntimeException("Simulation failed to find a valid " +
                    "Bennett parameter");
        }
        System.out.println("equilibration finished.");
        
        sim.integratorSim.getMoveManager().setEquilibrating(false);
        sim.setAccumulatorBlockSize(runBlockSize);
        sim.activityIntegrate.setMaxSteps(numSteps);
        sim.getController().actionPerformed();
        System.out.println("final reference optimal step frequency " + 
                sim.integratorSim.getStepFreq0() + " (actual: " + 
                sim.integratorSim.getActualStepFreq0() + ")");
        
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
        double ratio = sim.dsvo.getDataAsScalar();
        double error = sim.dsvo.getError();
        System.out.println("ratio average: "+ratio+", error: "+error);
        System.out.println("free energy difference: " + (-Math.log(ratio)) + 
                ", error: "+(error/ratio));
        System.out.println("target free energy: " + (AHarmonic-Math.log(ratio)));
        DataGroup allYourBase = 
            (DataGroup)sim.accumulators[0].getData(sim.dsvo.minDiffLocation());
        System.out.println("harmonic ratio average: " + 
                ((DataDoubleArray)allYourBase.getData(AccumulatorRatioAverage.StatType.RATIO.index)).getData()[1]
                 + " error: " + 
                ((DataDoubleArray)allYourBase.getData(AccumulatorRatioAverage.StatType.RATIO_ERROR.index)).getData()[1]);
        
        allYourBase = (DataGroup)sim.accumulators[1].getData(sim.accumulators[1].getNBennetPoints() -
                sim.dsvo.minDiffLocation()-1);
        System.out.println("target ratio average: " + 
                ((DataDoubleArray)allYourBase.getData(AccumulatorRatioAverage.StatType.RATIO.index)).getData()[1]
                 + " error: " + 
                ((DataDoubleArray)allYourBase.getData(AccumulatorRatioAverage.StatType.RATIO_ERROR.index)).getData()[1]);

        if(D==1) {
            double AHR = -(numMolecules-1)*Math.log(numMolecules/density-numMolecules)
                + SpecialFunctions.lnFactorial(numMolecules) ;
            System.out.println("Hard-rod free energy: "+AHR);
        }
        
        sim.sWriter.actionPerformed();
    }
    
    
    public void setComparedWV(int awv){
        compareMove.setComparedWV(awv);
        changeMove.setHarmonicWV(awv);
        meterBinA.setComparedWV(awv);
        meterBinB.setComparedWV(awv);
    }
    public static class SimOverlapSingleWaveVector3DParam extends ParameterBase {
        public int numAtoms = 32;
        public double density = 1.3;
        public int D = 3;
        public double harmonicFudge = 1.0;
        public String filename = "normal_modes_LJ_3D_32";
        public double temperature = 1.0;
        public int comparedWV = 7;
        
        public int numSteps = 40000000;
        public int runBlockSize = 100000;
        public int subBlockSize = 1000;    //# of steps in subintegrator per integrator step

        public int eqNumSteps = 4000000;  
        public int eqBlockSize = 10000;
        
        public int bennettNumSteps = 4000000;
        public int benBlockSize = 10000;

    }
    
    
    public static class MeterPotentialEnergyDifference extends DataSourceScalar {

        public MeterPotentialEnergyDifference(DataSourceScalar meter, double offset) {
            super("energy", Energy.DIMENSION);
            this.meter = meter;
            this.offset = offset;
        }
        
        public double getDataAsScalar() {
            return meter.getDataAsScalar() - offset;
        }

        private static final long serialVersionUID = 1L;
        protected final DataSourceScalar meter;
        protected final double offset;
    }
}
