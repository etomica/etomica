package etomica.models.oneDHardRods;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;

import etomica.action.activity.ActivityIntegrate;
import etomica.api.IAtomTypeLeaf;
import etomica.api.IBox;
import etomica.box.Box;
import etomica.data.AccumulatorAverageFixed;
import etomica.data.AccumulatorRatioAverage;
import etomica.data.DataPump;
import etomica.data.IEtomicaDataSource;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.data.types.DataDoubleArray;
import etomica.data.types.DataGroup;
import etomica.exception.ConfigurationOverlapException;
import etomica.integrator.IntegratorMC;
import etomica.lattice.crystal.BasisMonatomic;
import etomica.lattice.crystal.Primitive;
import etomica.lattice.crystal.PrimitiveCubic;
import etomica.math.SpecialFunctions;
import etomica.nbr.list.PotentialMasterList;
import etomica.normalmode.CoordinateDefinitionLeaf;
import etomica.normalmode.NormalModes1DHR;
import etomica.normalmode.P2XOrder;
import etomica.normalmode.WaveVectorFactory;
import etomica.potential.P2HardSphere;
import etomica.potential.Potential2;
import etomica.potential.Potential2HardSpherical;
import etomica.simulation.Simulation;
import etomica.space.Boundary;
import etomica.space.BoundaryRectangularPeriodic;
import etomica.space.Space;
import etomica.species.SpeciesSpheresMono;
import etomica.units.Null;
import etomica.util.ParameterBase;
import etomica.util.ReadParameters;
import etomica.virial.overlap.AccumulatorVirialOverlapSingleAverage;
import etomica.virial.overlap.DataSourceVirialOverlap;
import etomica.virial.overlap.IntegratorOverlap;

public class SimOverlapSingleWaveVector extends Simulation {
    private static final long serialVersionUID = 1L;
    private static final String APP_NAME = "SimOverlapAB";
    Primitive primitive;
    int[] nCells;
    NormalModes1DHR nm;
    double bennettParam;       //adjustable parameter - Bennett's parameter
    public IntegratorOverlap integratorSim; //integrator for the whole simulation
    public DataSourceVirialOverlap dsvo;
    public BasisMonatomic basis;
    ActivityIntegrate activityIntegrate;
    
    IntegratorMC[] integrators;
    protected int blockSize,eqBlockSize;
    
    public AccumulatorVirialOverlapSingleAverage[] accumulators;
    public DataPump[] accumulatorPumps;
    public IEtomicaDataSource[] meters;
    public IBox boxTarget, boxRef;
    public Boundary boundaryTarget, boundaryRef;
    MCMoveChangeMode changeMove;
    MCMoveCompareSingleMode compareMove;
    MeterPotentialEnergy meterAinB, meterAinA;
    MeterCompareSingleModeBrute meterBinA, meterBinB;
    MeterCompareTest meterTestBinA, meterTestBinB;
    MeterOverlap meterOverlapInA, meterOverlapInB;
    
    public SimOverlapSingleWaveVector(Space _space, int numAtoms, double density, double 
            temperature, String filename, double harmonicFudge, int awv){
        super(_space, true);
        
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
        PotentialMasterList potentialMasterTarget = new 
            PotentialMasterList(this, space);
        boxTarget = new Box(this, space);
        addBox(boxTarget);
        boxTarget.setNMolecules(species, numAtoms);
        
        Potential2 p2 = new P2HardSphere(space, 1.0, true);
        p2 = new P2XOrder(space, (Potential2HardSpherical)p2);
        p2.setBox(boxTarget);
        potentialMasterTarget.addPotential(p2, new IAtomTypeLeaf[]
                {species.getLeafType(), species.getLeafType()});
        
        primitive = new PrimitiveCubic(space, 1.0/density);
        boundaryTarget = new BoundaryRectangularPeriodic(space, getRandom(),
                numAtoms/density);
        nCells = new int[]{numAtoms};
        boxTarget.setBoundary(boundaryTarget);
        
        CoordinateDefinitionLeaf coordinateDefinitionTarget = new 
                CoordinateDefinitionLeaf(this, boxTarget, primitive, space);
        coordinateDefinitionTarget.initializeCoordinates(nCells);
        
        double neighborRange = 1.01/density;
        potentialMasterTarget.setRange(neighborRange);
        //find neighbors now.  Don't hook up NeighborListManager since the
        //  neighbors won't change
        potentialMasterTarget.getNeighborManager(boxTarget).reset();
        
        IntegratorMC integratorTarget = new IntegratorMC(potentialMasterTarget, 
                random, temperature);
        integrators[1] = integratorTarget;
        integratorTarget.setBox(boxTarget);
        
        nm = new NormalModes1DHR(space.D());
        nm.setHarmonicFudge(harmonicFudge);
        nm.setTemperature(temperature);
        
        WaveVectorFactory waveVectorFactoryTarget = nm.getWaveVectorFactory();
        waveVectorFactoryTarget.makeWaveVectors(boxTarget);
        
        changeMove = new MCMoveChangeMode(potentialMasterTarget, random);
        integratorTarget.getMoveManager().addMCMove(changeMove);
        changeMove.setWaveVectors(waveVectorFactoryTarget.getWaveVectors());
        changeMove.setWaveVectorCoefficients(
                waveVectorFactoryTarget.getCoefficients());
        changeMove.setEigenVectors(nm.getEigenvectors(boxTarget));
        changeMove.setCoordinateDefinition(coordinateDefinitionTarget);
        changeMove.setBox((IBox)boxTarget);
        changeMove.setStepSizeMin(0.001);
        changeMove.setStepSize(0.01);
        
        meterAinA = new MeterPotentialEnergy(potentialMasterTarget);
        meterAinA.setBox(boxTarget);
        
        meterBinA = new MeterCompareSingleModeBrute("meterBinA", potentialMasterTarget, 
                coordinateDefinitionTarget, boxTarget);
        meterBinA.setEigenVectors(nm.getEigenvectors(boxTarget));
        meterBinA.setOmegaSquared(nm.getOmegaSquared(boxTarget));
        meterBinA.setTemperature(temperature);
        meterBinA.setWaveVectorCoefficients(waveVectorFactoryTarget.getCoefficients());
        meterBinA.setWaveVectors(waveVectorFactoryTarget.getWaveVectors());
        
//        meterTestBinA = new MeterConvertTest("meterBinA", potentialMasterTarget, 
//                coordinateDefinitionTarget, boxTarget);
//        meterTestBinA.setEigenVectors(nm.getEigenvectors(boxTarget));
//        meterTestBinA.setOmegaSquared(nm.getOmegaSquared(boxTarget));
//        meterTestBinA.setTemperature(temperature);
//        meterTestBinA.setWaveVectorCoefficients(waveVectorFactoryTarget.getCoefficients());
//        meterTestBinA.setWaveVectors(waveVectorFactoryTarget.getWaveVectors());
//        AccumulatorAverageFixed binaSink = new AccumulatorAverageFixed();
//        DataPump pumpTestBinA = new DataPump(meterTestBinA, binaSink);
//        integratorTarget.addIntervalAction(pumpTestBinA);
//        integratorTarget.setActionInterval(pumpTestBinA, 1000000);
//        integratorTarget.setIntervalActionPriority(pumpTestBinA, 198);
        
        meterOverlapInA = new MeterOverlap("MeterOverlapInA", Null.DIMENSION, 
                meterAinA, meterBinA, temperature);
        meters[1] = meterOverlapInA;
        
        potentialMasterTarget.getNeighborManager(boxTarget).reset();
        
        
        
        
//REFERENCE        
        //Set up REFERENCE system - System B - 0 - hybrid system
        PotentialMasterList potentialMasterRef = new 
            PotentialMasterList(this, space);
        boxRef = new Box(this, space);
        addBox(boxRef);
        boxRef.setNMolecules(species, numAtoms);
//        accumulators[1] = new AccumulatorVirialOverlapSingleAverage(10, 11, true);
        
        p2 = new P2HardSphere(space, 1.0, true);
        p2 = new P2XOrder(space, (Potential2HardSpherical)p2);
        p2.setBox(boxRef);
        potentialMasterRef.addPotential(p2, new IAtomTypeLeaf[]
                {species.getLeafType(), species.getLeafType()});
        
        primitive = new PrimitiveCubic(space, 1.0/density);
        boundaryRef = new BoundaryRectangularPeriodic(space, getRandom(),
                numAtoms/density);
        nCells = new int[]{numAtoms};
        boxRef.setBoundary(boundaryRef);
        
        CoordinateDefinitionLeaf coordinateDefinitionRef = new 
                CoordinateDefinitionLeaf(this, boxRef, primitive, space);
        coordinateDefinitionRef.initializeCoordinates(nCells);
        
        neighborRange = 1.01/density;
        potentialMasterRef.setRange(neighborRange);
        //find neighbors now.  Don't hook up NeighborListManager since the
        //  neighbors won't change
        potentialMasterRef.getNeighborManager(boxRef).reset();
        
        IntegratorMC integratorRef = new IntegratorMC(potentialMasterRef, 
                random, temperature);
        integratorRef.setBox(boxRef);
        integrators[0] = integratorRef;
        
        nm = new NormalModes1DHR(space.D());
        nm.setHarmonicFudge(harmonicFudge);
        nm.setTemperature(temperature);
        
        WaveVectorFactory waveVectorFactoryRef = nm.getWaveVectorFactory();
        waveVectorFactoryRef.makeWaveVectors(boxRef);
        
        compareMove = new MCMoveCompareSingleMode(potentialMasterRef, 
                random);
        integratorRef.getMoveManager().addMCMove(compareMove);
        compareMove.setWaveVectors(waveVectorFactoryRef.getWaveVectors());
        compareMove.setWaveVectorCoefficients(waveVectorFactoryRef.getCoefficients());
        compareMove.setOmegaSquared(nm.getOmegaSquared(boxRef), 
                waveVectorFactoryRef.getCoefficients());
        compareMove.setEigenVectors(nm.getEigenvectors(boxRef));
        compareMove.setCoordinateDefinition(coordinateDefinitionRef);
        compareMove.setTemperature(temperature);
        compareMove.setBox((IBox)boxRef);
        compareMove.setStepSizeMin(0.001);
        compareMove.setStepSize(0.01);
        
        meterAinB = new MeterPotentialEnergy(potentialMasterRef);
        meterAinB.setBox(boxRef);
       
        meterBinB = new MeterCompareSingleModeBrute(potentialMasterRef,
                coordinateDefinitionRef, boxRef);
        meterBinB.setCoordinateDefinition(coordinateDefinitionRef);
        meterBinB.setEigenVectors(nm.getEigenvectors(boxRef));
        meterBinB.setOmegaSquared(nm.getOmegaSquared(boxRef));
        meterBinB.setTemperature(temperature);
        meterBinB.setWaveVectorCoefficients(waveVectorFactoryRef.getCoefficients());
        meterBinB.setWaveVectors(waveVectorFactoryRef.getWaveVectors());
        integratorRef.setMeterPotentialEnergy(meterBinB);
        
//        meterTestBinB = new MeterConvertTest(potentialMasterRef,
//                coordinateDefinitionRef, boxRef);
//        meterTestBinB.setCoordinateDefinition(coordinateDefinitionRef);
//        meterTestBinB.setEigenVectors(nm.getEigenvectors(boxRef));
//        meterTestBinB.setOmegaSquared(nm.getOmegaSquared(boxRef));
//        meterTestBinB.setTemperature(temperature);
//        meterTestBinB.setWaveVectorCoefficients(waveVectorFactoryRef.getCoefficients());
//        meterTestBinB.setWaveVectors(waveVectorFactoryRef.getWaveVectors());
//        integratorRef.setMeterPotentialEnergy(meterTestBinB);
//        AccumulatorAverageFixed binbSink = new AccumulatorAverageFixed();
//        DataPump pumpTestBinB = new DataPump(meterTestBinB, binbSink);
//        integratorRef.addIntervalAction(pumpTestBinB);
//        integratorRef.setActionInterval(pumpTestBinB, 1000000);
//        integratorRef.setIntervalActionPriority(pumpTestBinB, 199);
        
        meterOverlapInB = new MeterOverlap("MeterOverlapInB", Null.DIMENSION, 
                meterBinB, meterAinB, temperature);
        meters[0] = meterOverlapInB;
        
        integratorRef.setBox(boxRef);
        potentialMasterRef.getNeighborManager(boxRef).reset();
        
        
//JOINT
        //Set up the rest of the joint stuff
        setComparedWaveVector(awv);
        
        integratorSim = new IntegratorOverlap(random, new 
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
    
    public void initBennettParameter(String fileName, long initSteps) {
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
            activityIntegrate.setMaxSteps(initSteps/2);
            getController().actionPerformed();
            getController().reset();

            setAccumulator(new AccumulatorVirialOverlapSingleAverage(41,true),0);
            setAccumulator(new AccumulatorVirialOverlapSingleAverage(41,false),1);
            setBennettParameter(1e40,40);
            activityIntegrate.setMaxSteps(initSteps);
            getController().actionPerformed();
            getController().reset();

            int newMinDiffLoc = dsvo.minDiffLocation();
            bennettParam = accumulators[0].getBennetAverage(newMinDiffLoc)
                /accumulators[1].getBennetAverage(newMinDiffLoc);
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
        accumulators[iBox].setBlockSize(blockSize);
        if (accumulatorPumps[iBox] == null) {
            accumulatorPumps[iBox] = new DataPump(meters[iBox], newAccumulator);
            integrators[iBox].addIntervalAction(accumulatorPumps[iBox]);
            integrators[iBox].setActionInterval(accumulatorPumps[iBox], 
                    boxRef.getLeafList().getAtomCount()*2);
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
        blockSize = newBlockSize;
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
    public void equilibrate(String fileName, long initSteps) {
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
            setAccumulator(new AccumulatorVirialOverlapSingleAverage(1,true),0);
            setAccumulator(new AccumulatorVirialOverlapSingleAverage(1,false),1);
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
        SimOverlapSingleWaveVectorParam params = new SimOverlapSingleWaveVectorParam();
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
        
        long numSteps = params.numSteps;
        long blockSize = params.blockSize;
        long subBlockSize = params.subBlockSize;
        
        long numEqSteps = params.eqNumSteps;
        long eqBlockSize = params.eqBlockSize;

        long numBennettSteps = params.bennettNumSteps;
        long benBlockSize = params.benBlockSize;
        
        String refFileName = args.length > 0 ? filename+"_ref" : null;
        
        System.out.println("Running Nancy's 1DHR simulation");
        System.out.println(numMolecules+" atoms at density "+density);
        System.out.println("harmonic fudge: "+harmonicFudge);
        System.out.println("temperature: ");
        System.out.println("compared wave vector: " + comparedWV);
        System.out.println("Total steps: "+numSteps+" , split into blocks of "+blockSize);
        System.out.println(subBlockSize+" steps in subintegrator, per step in  main integrator");
        System.out.println(numEqSteps+" equilibration steps, split into blocks of "+ eqBlockSize);
        System.out.println(numBennettSteps+" Bennett-only steps, split into blocks of "+benBlockSize);
        //        System.out.println((numSteps/subBlockSize)+" total steps of " + subBlockSize);
        System.out.println("output data to "+filename);


        
        //instantiate simulations!
        SimOverlapSingleWaveVector sim = new SimOverlapSingleWaveVector(Space.getInstance(D), numMolecules,
                density, temperature, filename, harmonicFudge, comparedWV);
        System.out.println("instantiated");
        
        //Divide out all the steps, so that the subpieces have the proper # of steps
        numSteps /= (int)subBlockSize;
        numEqSteps /= (int)subBlockSize;
        numBennettSteps /= (int)subBlockSize;
        
        //start simulation & equilibrate
        sim.integratorSim.getMoveManager().setEquilibrating(true);
        sim.integratorSim.setNumSubSteps((int)subBlockSize);
        
        sim.setAccumulatorBlockSize((int)benBlockSize);
        sim.initBennettParameter(filename, numBennettSteps);
        if(Double.isNaN(sim.bennettParam) || sim.bennettParam == 0 || 
                Double.isInfinite(sim.bennettParam)){
//            System.out.println("bennet info:  " + )
            throw new RuntimeException("Simulation failed to find a valid " +
                    "Bennett parameter");
        }
        sim.setAccumulatorBlockSize((int)eqBlockSize);
        sim.equilibrate(refFileName, numEqSteps);
        if(Double.isNaN(sim.bennettParam) || sim.bennettParam == 0 || 
                Double.isInfinite(sim.bennettParam)){
//            System.out.println("bennet info:  " + )
            throw new RuntimeException("Simulation failed to find a valid " +
                    "Bennett parameter");
        }
        System.out.println("equilibration finished.");
        
        sim.integratorSim.getMoveManager().setEquilibrating(false);
        sim.activityIntegrate.setMaxSteps(numSteps);
        sim.getController().actionPerformed();
        System.out.println("final reference optimal step frequency " + 
                sim.integratorSim.getStepFreq0() + " (actual: " + 
                sim.integratorSim.getActualStepFreq0() + ")");
        
        double[][] omega2 = sim.nm.getOmegaSquared(sim.boxTarget); 
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
        
        
    }
      
    
    public void setComparedWaveVector(int awv){
        compareMove.setComparedWaveVector(awv);
        meterBinA.setComparedWV(awv);
        meterBinB.setComparedWV(awv);
    }
    public static class SimOverlapSingleWaveVectorParam extends ParameterBase {
        public int numAtoms = 32;
        public double density = 0.5;
        public int D = 1;
        public double harmonicFudge = 1.0;
        public String filename = "HR1D_";
        public double temperature = 1.0;
        public int comparedWV = 10;
        
        public long numSteps = 40000000;
        public long blockSize = 100000;
        public long subBlockSize = 1000;    //# of steps in subintegrator per integrator step

        public long eqNumSteps = 4000000;  
        public long eqBlockSize = 10000;
        
        public long bennettNumSteps = 2000000;
        public long benBlockSize = 10000;

    }
    
}
