package etomica.normalmode;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;

import etomica.action.activity.ActivityIntegrate;
import etomica.api.IAction;
import etomica.api.IAtomType;
import etomica.api.IBox;
import etomica.box.Box;
import etomica.data.AccumulatorAverage;
import etomica.data.AccumulatorAverageCollapsing;
import etomica.data.AccumulatorAverageFixed;
import etomica.data.AccumulatorHistogram;
import etomica.data.DataFork;
import etomica.data.DataLogger;
import etomica.data.DataPump;
import etomica.data.DataTableWriter;
import etomica.data.IEtomicaDataSource;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.data.meter.MeterPressure;
import etomica.data.types.DataDoubleArray;
import etomica.data.types.DataGroup;
import etomica.integrator.IntegratorBox;
import etomica.integrator.IntegratorMC;
import etomica.integrator.mcmove.MCMoveStepTracker;
import etomica.lattice.crystal.Basis;
import etomica.lattice.crystal.BasisCubicFcc;
import etomica.lattice.crystal.BasisMonatomic;
import etomica.lattice.crystal.Primitive;
import etomica.lattice.crystal.PrimitiveCubic;
import etomica.math.SpecialFunctions;
import etomica.potential.P2SoftSphere;
import etomica.potential.P2SoftSphericalTruncatedShifted;
import etomica.potential.Potential2SoftSpherical;
import etomica.potential.PotentialMasterMonatomic;
import etomica.simulation.Simulation;
import etomica.space.Boundary;
import etomica.space.BoundaryRectangularPeriodic;
import etomica.space.Space;
import etomica.species.SpeciesSpheresMono;
import etomica.util.DoubleRange;
import etomica.util.HistogramSimple;
import etomica.util.ParameterBase;
import etomica.util.ReadParameters;
import etomica.virial.overlap.AccumulatorVirialOverlapSingleAverage;
import etomica.virial.overlap.DataSourceVirialOverlap;
import etomica.virial.overlap.IntegratorOverlap;

/**
 * 
 * Simulation to run sampling with the soft sphere potential, but measuring
 * the harmonic potential based on normal mode data from a previous simulation.
 * 
 *  	Direct Sampling and Bennett's Overlap for Soft Sphere FCC structure
 *  
 * 
 * @author Andrew Schultz & Tai Tan
 */
public class SimDSBenOverlapSoftSphereFCC extends Simulation {

    public SimDSBenOverlapSoftSphereFCC(Space _space, int numAtoms, double density, double temperature, String filename, double harmonicFudge, int exponent) {
        super(_space, true);

        potentialMasterTarget = new PotentialMasterMonatomic(this);
        integrators = new IntegratorBox[2];
        accumulatorPumps = new DataPump[2];
        meters = new IEtomicaDataSource[2];
        accumulators = new AccumulatorVirialOverlapSingleAverage[2];

        SpeciesSpheresMono species = new SpeciesSpheresMono(this, space);
        getSpeciesManager().addSpecies(species);

        // TARGET
        boxTarget = new Box(space);
        addBox(boxTarget);
        boxTarget.setNMolecules(species, numAtoms);

        integratorTarget = new IntegratorMC(potentialMasterTarget, getRandom(), temperature);
        MCMoveAtomCoupled atomMove = new MCMoveAtomCoupled(potentialMasterTarget, getRandom(), space);
        atomMove.setStepSize(0.1);
        atomMove.setStepSizeMax(0.5);
        integratorTarget.getMoveManager().addMCMove(atomMove);
        ((MCMoveStepTracker)atomMove.getTracker()).setNoisyAdjustment(true);

        integrators[1] = integratorTarget;

        if (space.D() == 1) {
            primitive = new PrimitiveCubic(space, 1.0/density);
            boundaryTarget = new BoundaryRectangularPeriodic(space, numAtoms/density);
            nCells = new int[]{numAtoms};
            basis = new BasisMonatomic(space);
        } else {
            double L = Math.pow(4.0/density, 1.0/3.0);
            primitive = new PrimitiveCubic(space, L);
            int n = (int)Math.round(Math.pow(numAtoms/4, 1.0/3.0));
            nCells = new int[]{n,n,n};
            boundaryTarget = new BoundaryRectangularPeriodic(space, n * L);
            basis = new BasisCubicFcc();
        }
        boxTarget.setBoundary(boundaryTarget);

        CoordinateDefinitionLeaf coordinateDefinitionTarget = new CoordinateDefinitionLeaf(this, boxTarget, primitive, basis, space);
        coordinateDefinitionTarget.initializeCoordinates(nCells);

        Potential2SoftSpherical potential = new P2SoftSphere(space, 1.0, 1.0, exponent);
        double truncationRadius = boundaryTarget.getDimensions().x(0) * 0.495;
        P2SoftSphericalTruncatedShifted pTruncated = new P2SoftSphericalTruncatedShifted(space, potential, truncationRadius);
        IAtomType sphereType = species.getLeafType();
        potentialMasterTarget.addPotential(pTruncated, new IAtomType[] { sphereType, sphereType });
        atomMove.setPotential(pTruncated);

        integratorTarget.setBox(boxTarget);

        /*
         *  1-body Potential to Constraint the atom from moving too far
         *  	away from its lattice-site
         *  
         */

        P1Constraint p1Constraint = new P1Constraint(space, primitive, boxTarget, coordinateDefinitionTarget);
        potentialMasterTarget.addPotential(p1Constraint, new IAtomType[] {sphereType});
        
        potentialMasterTarget.lrcMaster().setEnabled(false);
        MeterPotentialEnergy meterPE = new MeterPotentialEnergy(potentialMasterTarget);
        meterPE.setBox(boxTarget);
        latticeEnergy = meterPE.getDataAsScalar();
        
    
        // HARMONIC
        boundaryHarmonic = new BoundaryRectangularPeriodic(space);
        boxHarmonic = new Box(boundaryHarmonic, space);
        addBox(boxHarmonic);
        boxHarmonic.setNMolecules(species, numAtoms);

        IntegratorMC integratorHarmonic = new IntegratorMC(potentialMasterTarget, random, 1.0);

        move = new MCMoveHarmonic(getRandom());
        integratorHarmonic.getMoveManager().addMCMove(move);
        integrators[0] = integratorHarmonic;
        
        if (space.D() == 1) {
            boundaryHarmonic = new BoundaryRectangularPeriodic(space, numAtoms/density);
        } else {
            double L = Math.pow(4.0/density, 1.0/3.0);
            int n = (int)Math.round(Math.pow(numAtoms/4, 1.0/3.0));
            boundaryHarmonic = new BoundaryRectangularPeriodic(space, n * L);
        }
        boxHarmonic.setBoundary(boundaryHarmonic);

        CoordinateDefinitionLeaf coordinateDefinitionHarmonic = new CoordinateDefinitionLeaf(this, boxHarmonic, primitive, basis, space);
        coordinateDefinitionHarmonic.initializeCoordinates(nCells);
        
        normalModes = new NormalModesFromFile(filename, space.D());
        normalModes.setHarmonicFudge(harmonicFudge);
        
        /*
         * nuke this line if it is overlap between DB and harmonic
         */
       normalModes.setTemperature(temperature);
        
        WaveVectorFactory waveVectorFactory = normalModes.getWaveVectorFactory();
        waveVectorFactory.makeWaveVectors(boxHarmonic);
        move.setOmegaSquared(normalModes.getOmegaSquared(boxHarmonic), waveVectorFactory.getCoefficients());
        move.setEigenVectors(normalModes.getEigenvectors(boxHarmonic));
        move.setWaveVectors(waveVectorFactory.getWaveVectors());
        move.setWaveVectorCoefficients(waveVectorFactory.getCoefficients());
        move.setCoordinateDefinition(coordinateDefinitionHarmonic);
        move.setTemperature(temperature);
        
        move.setBox(boxHarmonic);
        
        integratorHarmonic.setBox(boxHarmonic);

        // OVERLAP
        integratorOverlap = new IntegratorOverlap(random, new IntegratorBox[]{integratorHarmonic, integratorTarget});
        //integratorOverlap.setAdjustStepFreq(false);
        //integratorOverlap.setStepFreq0(0.02); //to make the number of sampling the same!! 7/29/08
        
        meterHarmonicEnergy = new MeterHarmonicEnergy(coordinateDefinitionTarget, normalModes);
        meterHarmonicEnergy.setBox(boxTarget);
        
        // target ---> harmonic
        MeterBoltzmannTarget meterTarget = new MeterBoltzmannTarget(integratorTarget, meterHarmonicEnergy);
        meterTarget.setLatticeEnergy(latticeEnergy);
        meters[1] = meterTarget;
        setAccumulator(new AccumulatorVirialOverlapSingleAverage(10, 11, false), 1);
        
        // harmonic ---> target
        MeterBoltzmannHarmonic meterHarmonic = new MeterBoltzmannHarmonic(move, potentialMasterTarget);
        meterHarmonic.setTemperature(temperature);
        meterHarmonic.setLatticeEnergy(latticeEnergy);
        meters[0] = meterHarmonic;
        setAccumulator(new AccumulatorVirialOverlapSingleAverage(10, 11, true), 0);
        
        setRefPref(1.0, 30);
     
        activityIntegrate = new ActivityIntegrate(integratorOverlap);
        
        getController().addAction(activityIntegrate);
    }

    public void setRefPref(double refPrefCenter, double span) {
        refPref = refPrefCenter;
        accumulators[0].setBennetParam(refPrefCenter,span);
        accumulators[1].setBennetParam(refPrefCenter,span);
        // needed for Bennett sampling
//        if (accumulators[0].getNBennetPoints() == 1) {
//            ((MeterBoltzmannHarmonic)meters[0]).refPref = refPrefCenter;
//            ((MeterBoltzmannTarget)meters[1]).refPref = refPrefCenter;
//        }
    }

    public void setAccumulator(AccumulatorVirialOverlapSingleAverage newAccumulator, int iBox) {
        accumulators[iBox] = newAccumulator;
        newAccumulator.setBlockSize(1);
        if (accumulatorPumps[iBox] == null) {
            accumulatorPumps[iBox] = new DataPump(meters[iBox],newAccumulator);
            integrators[iBox].addIntervalAction(accumulatorPumps[iBox]);
            if (iBox == 1) {
            	if (boxTarget.getMoleculeList().getMoleculeCount()==32){
            		
            		integrators[iBox].setActionInterval(accumulatorPumps[iBox], 500);
            	
            	} else if (boxTarget.getMoleculeList().getMoleculeCount()==108){
                
            		integrators[iBox].setActionInterval(accumulatorPumps[iBox], 1000);
            	} else 
            		
            		integrators[iBox].setActionInterval(accumulatorPumps[iBox], boxTarget.getMoleculeList().getMoleculeCount());
            }
        }
        else {
            accumulatorPumps[iBox].setDataSink(newAccumulator);
        }
        if (integratorOverlap != null && accumulators[0] != null && accumulators[1] != null) {
            dsvo = new DataSourceVirialOverlap(accumulators[0],accumulators[1]);
            integratorOverlap.setDSVO(dsvo);
        }
    }
    
    public void setRefPref(double newRefPref) {
        System.out.println("setting ref pref (explicitly) to "+newRefPref);
        setAccumulator(new AccumulatorVirialOverlapSingleAverage(1, 1,true),0);
        setAccumulator(new AccumulatorVirialOverlapSingleAverage(1,false),1);
        setRefPref(newRefPref,1);
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
                System.out.println("setting ref pref (from file) to "+refPref);
                setAccumulator(new AccumulatorVirialOverlapSingleAverage(1,1,true),0);
                setAccumulator(new AccumulatorVirialOverlapSingleAverage(1,false),1);
                setRefPref(refPref,1);
            }
            catch (IOException e) {
                // file not there, which is ok.
            }
        }
        
        if (refPref == -1) {
            // equilibrate off the lattice to avoid anomolous contributions
            activityIntegrate.setMaxSteps(initSteps/2);
            getController().actionPerformed();
            getController().reset();
            System.out.println("target equilibration finished");

            setAccumulator(new AccumulatorVirialOverlapSingleAverage(1,41,true),0);
            setAccumulator(new AccumulatorVirialOverlapSingleAverage(41,false),1);
            setRefPref(1.0,200);
            activityIntegrate.setMaxSteps(initSteps);
            getController().actionPerformed();
            getController().reset();

            int newMinDiffLoc = dsvo.minDiffLocation();
            refPref = accumulators[0].getBennetAverage(newMinDiffLoc)
                /accumulators[1].getBennetAverage(newMinDiffLoc);
            if (Double.isNaN(refPref) || refPref == 0 || Double.isInfinite(refPref)) {
                throw new RuntimeException("Simulation failed to find a valid ref pref");
            }
            System.out.println("setting ref pref to "+refPref);
            
            setAccumulator(new AccumulatorVirialOverlapSingleAverage(1,11,true),0);
            setAccumulator(new AccumulatorVirialOverlapSingleAverage(11,false),1);
            setRefPref(refPref,5);

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

        for (int i=0; i<2; i++) {
            if (integrators[i] instanceof IntegratorMC) ((IntegratorMC)integrators[i]).getMoveManager().setEquilibrating(true);
        }
        getController().actionPerformed();
        getController().reset();
        for (int i=0; i<2; i++) {
            if (integrators[i] instanceof IntegratorMC) ((IntegratorMC)integrators[i]).getMoveManager().setEquilibrating(false);
        }

        if (refPref == -1) {
            int newMinDiffLoc = dsvo.minDiffLocation();
            refPref = accumulators[0].getBennetAverage(newMinDiffLoc)
                /accumulators[1].getBennetAverage(newMinDiffLoc);
            System.out.println("setting ref pref to "+refPref+" ("+newMinDiffLoc+")");
            setAccumulator(new AccumulatorVirialOverlapSingleAverage(1,1,true),0);
            
            System.out.println("block size (equilibrate) "+accumulators[0].getBlockSize());
            
            setAccumulator(new AccumulatorVirialOverlapSingleAverage(1,false),1);
            setRefPref(refPref,1);
            
            if (fileName != null) {
                try {
                    FileWriter fileWriter = new FileWriter(fileName);
                    BufferedWriter bufWriter = new BufferedWriter(fileWriter);
                    bufWriter.write(String.valueOf(refPref)+"\n");
                    bufWriter.close();
                    fileWriter.close();
                }
                catch (IOException e) {
                    throw new RuntimeException("couldn't write to refpref file");
                }
            }
        }
        else {
            dsvo.reset();
        }
    }

    /**
     * @param args filename containing simulation parameters
     * @see SimDSBenOverlapSoftSphereFCC.SimOverlapParam
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
        int numMolecules = params.numMolecules;
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
        final SimDSBenOverlapSoftSphereFCC sim = new SimDSBenOverlapSoftSphereFCC(Space.getInstance(D), numMolecules, density, temperature, filename, harmonicFudge, exponentN);
        
        //start simulation
        sim.integratorOverlap.setNumSubSteps(1000);   
        numSteps /= 1000;

//        StopWatcher stopWatcher = new StopWatcher("stop", sim, filename);
//        IntervalActionAdapter iaa = new IntervalActionAdapter(stopWatcher);
//        iaa.setActionInterval(100);
//        sim.integratorOverlap.addListener(iaa);

        //sim.integratorOverlap.setAdjustStepFreq(false);
        //sim.integratorOverlap.setStepFreq0(0.02);  //11/5/2008

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

        System.out.println("final reference optimal step frequency "+sim.integratorOverlap.getStepFreq0()
        		+" (actual: "+sim.integratorOverlap.getActualStepFreq0()+")");

        IEtomicaDataSource[] workMeters = new IEtomicaDataSource[2];
        IEtomicaDataSource[] workBennets = new IEtomicaDataSource[2];
        IEtomicaDataSource[] boltzmannDirectSampling = new IEtomicaDataSource[2];
      
      /*
       *	 Harmonic Sampling
       */
        
        //start of Harmonic
        /*
         * Direct Sampling
         * Boltzmann Factor difference average (Harmonic ---> Target)  
         */
        MeterDirectSamplingHarmonic meterDirectSamplingHarmonic = new MeterDirectSamplingHarmonic(sim.move, sim.potentialMasterTarget);
        meterDirectSamplingHarmonic.setTemperature(temperature);
        meterDirectSamplingHarmonic.setLatticeEnergy(sim.latticeEnergy);
        
        boltzmannDirectSampling[0] = meterDirectSamplingHarmonic;
        
        final AccumulatorAverageFixed dataAverageBoltzmannHarmonic = new AccumulatorAverageFixed(1);
        DataPump pumpBoltzmannHarmonic = new DataPump(boltzmannDirectSampling[0], dataAverageBoltzmannHarmonic);
        
        sim.integrators[0].addIntervalAction(pumpBoltzmannHarmonic);
        sim.integrators[0].setActionInterval(pumpBoltzmannHarmonic, 1);
        
        // Work Harmonic --> Target
        MeterWorkHarmonicPhaseSpace meterWorkHarmonic = new MeterWorkHarmonicPhaseSpace(sim.move, sim.potentialMasterTarget);
        meterWorkHarmonic.setTemperature(temperature);
        meterWorkHarmonic.setLatticeEnergy(sim.latticeEnergy);
        workMeters[0] = meterWorkHarmonic;
        
        DataFork dataForkHarmonic = new DataFork();
        DataPump pumpHarmonic = new DataPump(workMeters[0], dataForkHarmonic);
        
        final AccumulatorAverageFixed dataAverageHarmonic = new AccumulatorAverageFixed(1);
        dataForkHarmonic.addDataSink(dataAverageHarmonic);
        sim.integrators[0].addIntervalAction(pumpHarmonic);
        sim.integrators[0].setActionInterval(pumpHarmonic, 1);
       
        //Histogram Work Harmonic ---> Target
        final AccumulatorHistogram histogramHarmonicTarget = new AccumulatorHistogram(new HistogramSimple(2500, new DoubleRange(-50, 200)));
        dataForkHarmonic.addDataSink(histogramHarmonicTarget);
        
        // end of Harmonic
        
        
     /*
      * 	Target Sampling
      */
        
        // start of Target
        /*
         * Direct Sampling
         * Boltzmann Factor difference average (Target ---> Harmonic)
         */
        MeterDirectSamplingTarget meterDirectSamplingTarget = new MeterDirectSamplingTarget(sim.integratorTarget, sim.meterHarmonicEnergy);
        meterDirectSamplingTarget.setLatticeEnergy(sim.latticeEnergy);
        boltzmannDirectSampling[1] = meterDirectSamplingTarget;
        
        final AccumulatorAverageFixed dataAverageBoltzmannTarget = new AccumulatorAverageFixed();
        DataPump pumpBoltzmannTarget = new DataPump(boltzmannDirectSampling[1], dataAverageBoltzmannTarget);
        
        sim.integrators[1].addIntervalAction(pumpBoltzmannTarget);
        sim.integrators[1].setActionInterval(pumpBoltzmannTarget, 1);
        
        // Work Target --> Harmonic
        MeterWorkTargetPhaseSpace meterWorkTarget =  new MeterWorkTargetPhaseSpace(sim.integratorTarget, sim.meterHarmonicEnergy);
        meterWorkTarget.setLatticeEnergy(sim.latticeEnergy);
        workMeters[1] = meterWorkTarget;
        
        DataFork dataForkTarget = new DataFork();
        DataPump pumpTarget = new DataPump(workMeters[1], dataForkTarget);
        
        final AccumulatorAverageFixed dataAverageTarget = new AccumulatorAverageFixed(1);
        dataForkTarget.addDataSink(dataAverageTarget);
        sim.integrators[1].addIntervalAction(pumpTarget);
        sim.integrators[1].setActionInterval(pumpTarget, 1);
        
        //Histogram Work Target ---> Harmonic
        final AccumulatorHistogram histogramTargetHarmonic = new AccumulatorHistogram(new HistogramSimple(2500, new DoubleRange(-50, 200)));
        dataForkTarget.addDataSink(histogramTargetHarmonic);
        
        //Pressure (Integrator Target)
        MeterPressure meterPressureTarget = new MeterPressure(Space.getInstance(D));
        meterPressureTarget.setIntegrator(sim.integrators[1]);
        
        final AccumulatorAverage pressureTargetAverage = new AccumulatorAverageCollapsing();
        DataPump pumpPressureTarget = new DataPump(meterPressureTarget, pressureTargetAverage);
        sim.integrators[1].addIntervalAction(pumpPressureTarget);
        sim.integrators[1].setActionInterval(pumpPressureTarget, 100);
        // end of Target
    
        
   /*
    * Perturbation into the overlap Region
    */
      
        // Harmonic ---> Bennet's Overlap
         final MeterWorkHarmonicBennet meterWorkHarmonicBennet = new MeterWorkHarmonicBennet(sim.move, sim.potentialMasterTarget, sim.refPref);
         meterWorkHarmonicBennet.setTemperature(temperature);
         meterWorkHarmonicBennet.setLatticeEnergy(sim.latticeEnergy);
         workBennets[0] = meterWorkHarmonicBennet;
          
         DataFork dataForkHarmonicBennet = new DataFork();
         DataPump pumpHarmonicBennet = new DataPump(workBennets[0], dataForkHarmonicBennet);
          
         final AccumulatorAverageFixed dataAverageHarmonicBennet = new AccumulatorAverageFixed(1);
          
         dataForkHarmonicBennet.addDataSink(dataAverageHarmonicBennet);
         sim.integrators[0].addIntervalAction(pumpHarmonicBennet);
         sim.integrators[0].setActionInterval(pumpHarmonicBennet, 1);
        
         //Histogram Harmonic--> Bennett's
         final AccumulatorHistogram histogramHarmonicBenn = new AccumulatorHistogram(new HistogramSimple(2500, new DoubleRange(-50,200)));
         dataForkHarmonicBennet.addDataSink(histogramHarmonicBenn);
        
      /*
       *  Target ---> Bennet's Overlap
       */
         final MeterWorkTargetBennet meterWorkTargetBennet = new MeterWorkTargetBennet(sim.integratorTarget, sim.meterHarmonicEnergy,sim.refPref);
         meterWorkTargetBennet.setLatticeEnergy(sim.latticeEnergy);
         workBennets[1] = meterWorkTargetBennet;
        
         DataFork dataForkTargetBennet = new DataFork();
         DataPump pumpTargetBennet = new DataPump(workBennets[1], dataForkTargetBennet);
        
         final AccumulatorAverageFixed dataAverageTargetBennet = new AccumulatorAverageFixed();
        
         dataForkTargetBennet.addDataSink(dataAverageTargetBennet);
         sim.integrators[1].addIntervalAction(pumpTargetBennet);
         sim.integrators[1].setActionInterval(pumpTargetBennet, numMolecules*2);
        
         //Histogram Target--> Bennett's
         final AccumulatorHistogram histogramTargetBenn = new AccumulatorHistogram(new HistogramSimple(2500,new DoubleRange(-50,200)));
         dataForkTargetBennet.addDataSink(histogramTargetBenn);


        
        FileWriter fileWriter, fileWriterBennet, 
        		   fileWriterHarmonic, fileWriterTarget;
               
        try{
        	fileWriter = new FileWriter(filename+"_DSFE");
        	fileWriterBennet = new FileWriter(filename+"_Ben");
        	fileWriterHarmonic = new FileWriter(filename+ "_Harm");
        	fileWriterTarget = new FileWriter(filename+ "_Targ");
        }catch(IOException e){
        	fileWriter = null;
        	fileWriterBennet = null;
        	fileWriterHarmonic = null;
        	fileWriterTarget = null;
        }
        
        final double temp = temperature;
        final int dimension = D;
        final String outFileName = filename;
        final FileWriter fileWriterDSFE = fileWriter;
        final FileWriter fileWriterBen = fileWriterBennet;
        final FileWriter fileWriterHarm = fileWriterHarmonic;
        final FileWriter fileWriterTarg = fileWriterTarget;
        
        IAction outputActionOverlap = new IAction(){
        	public void actionPerformed(){
        		long idStep = sim.integratorOverlap.getStepCount();
     
		        /*
		         * Histogram
		         */
		        //Target ---> Harmonic
				DataLogger dataLogger1 = new DataLogger();
				DataTableWriter dataTableWriter1 = new DataTableWriter();
				dataLogger1.setFileName(outFileName + "_hist_TargHarm");
				dataLogger1.setDataSink(dataTableWriter1);
				dataTableWriter1.setIncludeHeader(false);
				dataLogger1.putDataInfo(histogramTargetHarmonic.getDataInfo());
				
				dataLogger1.setWriteInterval(1);
				dataLogger1.setAppending(false); //overwrite data
				dataLogger1.putData(histogramTargetHarmonic.getData());
				
				dataLogger1.closeFile();
		       
				//Harmonic ---> Target
				DataLogger dataLogger2 = new DataLogger();
				DataTableWriter dataTableWriter2 = new DataTableWriter();
				dataLogger2.setFileName(outFileName + "_hist_HarmTarg");
		        dataLogger2.setDataSink(dataTableWriter2);
				dataTableWriter2.setIncludeHeader(false);
				dataLogger2.putDataInfo(histogramHarmonicTarget.getDataInfo());
		        
				dataLogger2.setWriteInterval(1);
				dataLogger2.setAppending(false); //overwrite data
		        dataLogger2.putData(histogramHarmonicTarget.getData());

		        dataLogger2.closeFile();
		        
		        
				//Harmonic ---> Bennett's 
				DataLogger dataLogger3 = new DataLogger();
				DataTableWriter dataTableWriter3 = new DataTableWriter();
				dataLogger3.setFileName(outFileName + "_hist_HarmBenn");
				dataLogger3.setDataSink(dataTableWriter3);
				dataTableWriter3.setIncludeHeader(false);
				dataLogger3.putDataInfo(histogramHarmonicBenn.getDataInfo());
				
				dataLogger3.setWriteInterval(1);
				dataLogger3.setAppending(false); //overwrite data
				dataLogger3.putData(histogramHarmonicBenn.getData());
				dataLogger3.closeFile();
			
		        
		        //Target ---> Bennet's
				DataLogger dataLogger4 = new DataLogger();
				DataTableWriter dataTableWriter4 = new DataTableWriter();
				dataLogger4.setFileName(outFileName + "_hist_TargBenn");
				dataLogger4.setDataSink(dataTableWriter4);
				dataTableWriter4.setIncludeHeader(false);
				dataLogger4.putDataInfo(histogramTargetBenn.getDataInfo());
				
				
				dataLogger4.setWriteInterval(1);
				dataLogger4.setAppending(false); //overwrite data
				dataLogger4.putData(histogramTargetBenn.getData());
				dataLogger4.closeFile();
		        
		        /*
		         * 
		         */

		        double[][] omega2 = sim.normalModes.getOmegaSquared(sim.boxTarget);
		        double[] coeffs = sim.normalModes.getWaveVectorFactory().getCoefficients();
		        double AHarmonic = 0;
		        for(int i=0; i<omega2.length; i++) {
		            for(int j=0; j<omega2[0].length; j++) {
		                if (!Double.isInfinite(omega2[i][j])) {
		                    AHarmonic += coeffs[i]*Math.log(omega2[i][j]*coeffs[i]/(temp*Math.PI));
		                }
		            }           
		            
		        }
		
		        int totalCells = 1;
		        for (int i=0; i<dimension; i++) {
		            totalCells *= sim.nCells[i];
		        }
		        int basisSize = sim.basis.getScaledCoordinates().length;
		        double fac = 1;
		        if (totalCells % 2 == 0) {
		            fac = Math.pow(2,dimension);
		        }
		        AHarmonic -= Math.log(Math.pow(2.0, basisSize*dimension*(totalCells - fac)/2.0) / Math.pow(totalCells,0.5*dimension));
		  
		        /*
		         * ratio = Q_target/Q_reference
		         * delta_FE = FE_target- FE_reference = - ln (ratio)
		         * 
		         * FE_target = FE_reference + deltaFE
		         * 
		         * free energy difference = temp * (FE_reference + FE_target)
		         * 
		         *  note: Aharmonic has not been timed temperature yet
		         *        that is why FE_reference * temp
		         */
		        double ratio = sim.dsvo.getDataAsScalar();
		        double error = sim.dsvo.getError();
		        double deltaFE = -Math.log(ratio);
		        double targetFE = temp*(AHarmonic + deltaFE);
		        double errorFE = temp*(error/ratio);
		        
		        /*
		         * Direct Sampling 
		         */
		        double boltzmannHarmonicAverage = dataAverageBoltzmannHarmonic.getData().getValue(AccumulatorAverage.StatType.AVERAGE.index);
		        double boltzmannHarmonicError = dataAverageBoltzmannHarmonic.getData().getValue(AccumulatorAverage.StatType.ERROR.index);
		        
		        double boltzmannTargetAverage = dataAverageBoltzmannTarget.getData().getValue(AccumulatorAverage.StatType.AVERAGE.index);
		        double boltzmannTargetError = dataAverageBoltzmannTarget.getData().getValue(AccumulatorAverage.StatType.ERROR.index);
		        
		        double deltaFE_DSHarmonic = -Math.log(boltzmannHarmonicAverage);
		        double deltaFE_DSTarget = -Math.log(boltzmannTargetAverage);
		        
		        double wHarmonic = dataAverageHarmonic.getData().getValue(AccumulatorAverage.StatType.AVERAGE.index);
		        double wTarget = dataAverageTarget.getData().getValue(AccumulatorAverage.StatType.AVERAGE.index);
		
		        double eHarmonic = dataAverageHarmonic.getData().getValue(AccumulatorAverage.StatType.ERROR.index);
		        double eTarget = dataAverageTarget.getData().getValue(AccumulatorAverage.StatType.ERROR.index);
		        
		        /*
		         * s_A = (u_B - u_A) - (FE_B - FE_A) 
		         * s_B = (u_A - u_B) + (FE_B - FE_A)
		         * 
		         * 	[  deltaFE = FE_B - FE_A  ]
		         *  
		         *  
		         *  A: target system
		         *  B: reference system
		         */
		        
		        double sHarmonic = wHarmonic - deltaFE;
		        double sTarget =   wTarget + deltaFE;
		        
		        double er_sHarmonic = Math.sqrt(eHarmonic*eHarmonic + (error/ratio)*(error/ratio));
		        double er_sTarget = Math.sqrt(eTarget*eTarget + (error/ratio)*(error/ratio));
		        
		        
		        /*
		         *     Q1/Q0 = ( Q_overlap/Q0 ) / ( Q_overlap/Q1 )    
		         *         
		         *         
		         *    s_A = (u_B - u_A) - (FE_B - FE_A) 
		         * 
		         *   	[  deltaFE_harmonic = FE_B - FE_A  ]
		         *      [  deltaFE_target   = FE_B - FE_A  ] 
		         *  
		         *  
		         *  A: reference system, target system
		         *  B: overlap region
		         */
		         
		        double wHarmonicBennet = dataAverageHarmonicBennet.getData().getValue(AccumulatorAverage.StatType.AVERAGE.index);
		        double wTargetBennet = dataAverageTargetBennet.getData().getValue(AccumulatorAverage.StatType.AVERAGE.index);
		
		        double eHarmonicBennet = dataAverageHarmonicBennet.getData().getValue(AccumulatorAverage.StatType.ERROR.index);
		        double eTargetBennet = dataAverageTargetBennet.getData().getValue(AccumulatorAverage.StatType.ERROR.index);
		        
		        
		        DataGroup allYourBase0 = (DataGroup)sim.accumulators[0].getData(sim.dsvo.minDiffLocation());
		        double ratioHarmonicAverage = ((DataDoubleArray)allYourBase0.getData(AccumulatorAverage.StatType.AVERAGE.index)).getData()[1];
		        double ratioHarmonicError = ((DataDoubleArray)allYourBase0.getData(AccumulatorAverage.StatType.ERROR.index)).getData()[1];
		        double deltaFEHarmonic = -Math.log(ratioHarmonicAverage);

		        DataGroup allYourBase1 = (DataGroup)sim.accumulators[1].getData(sim.dsvo.minDiffLocation());
		        double ratioTargetAverage = ((DataDoubleArray)allYourBase1.getData(AccumulatorAverage.StatType.AVERAGE.index)).getData()[1];
		        double ratioTargetError = ((DataDoubleArray)allYourBase1.getData(AccumulatorAverage.StatType.ERROR.index)).getData()[1];
		        double deltaFETarget = -Math.log(ratioTargetAverage);
		        
		        double sHarmonicBennet = wHarmonicBennet - deltaFEHarmonic;
		        double sTargetBennet = wTargetBennet - deltaFETarget;
		        
		        /*
		         *  Reversed direction (Bennett's to [harmonic or target])
		         *  reweighted work determined
		         */
		        double reweightedWorkBennHarm = meterWorkHarmonicBennet.getDataReweighted();
		        double reweightedWorkBennTarg = meterWorkTargetBennet.getDataReweighted();
		        
		        double sBennetHarmonic = -reweightedWorkBennHarm + deltaFEHarmonic;
		        double sBennetTarget = -reweightedWorkBennTarg + deltaFETarget;
		        
		        double pressureTarget = ((DataGroup)pressureTargetAverage.getData()).getValue(AccumulatorAverage.StatType.AVERAGE.index);
		        double pressureError = ((DataGroup)pressureTargetAverage.getData()).getValue(AccumulatorAverage.StatType.ERROR.index);
		        
		        try {
		        	
		        	fileWriterDSFE.write(idStep*1000 + " " + deltaFE + " " + deltaFE_DSHarmonic + " " + deltaFE_DSTarget 
		        							  + " " + targetFE + " "+ errorFE + " "+ wHarmonic + " " + wTarget
		        			                  + " " + sHarmonic+ " "+ sTarget + " "+ pressureTarget + " " + pressureError + "\n");
		        	
		        	fileWriterBen.write(idStep*1000 + " " + deltaFEHarmonic + " " + deltaFETarget + " " + wHarmonicBennet + " "+ wTargetBennet
		        			+ " " + reweightedWorkBennHarm + " " + reweightedWorkBennTarg                      
		        			+ " " + sHarmonicBennet + " " + sTargetBennet 
		        			+ " " + sBennetHarmonic + " " + sBennetTarget + "\n");
		        	
		        } catch (IOException e){
		        	
		        }
        	}
        }; 
        
        
        
        IAction outputActionHarmonic = new IAction(){
        	public void actionPerformed(){

        		long idStep = sim.integrators[0].getStepCount();
        		
				/*
				 * Histogram
				 */
				//Harmonic ---> Target
				DataLogger dataLogger1 = new DataLogger();
				DataTableWriter dataTableWriter1 = new DataTableWriter();
				dataLogger1.setFileName(outFileName + "_histDS_HarmTarg");
		        dataLogger1.setDataSink(dataTableWriter1);
				dataTableWriter1.setIncludeHeader(false);
				dataLogger1.putDataInfo(histogramHarmonicTarget.getDataInfo());
		        
				dataLogger1.setWriteInterval(1);
				dataLogger1.setAppending(false); //overwrite data
		        dataLogger1.putData(histogramHarmonicTarget.getData());

		        dataLogger1.closeFile();

		        double[][] omega2 = sim.normalModes.getOmegaSquared(sim.boxTarget);
		        double[] coeffs = sim.normalModes.getWaveVectorFactory().getCoefficients();
		        double AHarmonic = 0;
		        for(int i=0; i<omega2.length; i++) {
		            for(int j=0; j<omega2[0].length; j++) {
		                if (!Double.isInfinite(omega2[i][j])) {
		                    AHarmonic += coeffs[i]*Math.log(omega2[i][j]*coeffs[i]/(temp*Math.PI));
		                }
		            }           
		            
		        }
		
		        int totalCells = 1;
		        for (int i=0; i<dimension; i++) {
		            totalCells *= sim.nCells[i];
		        }
		        int basisSize = sim.basis.getScaledCoordinates().length;
		        double fac = 1;
		        if (totalCells % 2 == 0) {
		            fac = Math.pow(2,dimension);
		        }
		        AHarmonic -= Math.log(Math.pow(2.0, basisSize*dimension*(totalCells - fac)/2.0) / Math.pow(totalCells,0.5*dimension));
 
		        /*
		         * ratio = Q_target/Q_reference
		         * delta_FE = FE_target- FE_reference = - ln (ratio)
		         * 
		         * FE_target = FE_reference + deltaFE
		         * 
		         * free energy difference = temp * (FE_reference + FE_target)
		         * 
		         *  note: Aharmonic has not been timed temperature yet
		         *        that is why FE_reference * temp
		         */
		        double ratio = sim.dsvo.getDataAsScalar();
		        double deltaFE = -Math.log(ratio);
		        
		        double boltzmannHarmonicAverage = dataAverageBoltzmannHarmonic.getData().getValue(AccumulatorAverage.StatType.AVERAGE.index);
		        double deltaFE_DSHarmonic = -Math.log(boltzmannHarmonicAverage);
		        
		        double wHarmonic = dataAverageHarmonic.getData().getValue(AccumulatorAverage.StatType.AVERAGE.index);
		        double wTarget = dataAverageTarget.getData().getValue(AccumulatorAverage.StatType.AVERAGE.index);
		
		        /*
		         * s_A = (u_B - u_A) - (FE_B - FE_A) 
		         * s_B = (u_A - u_B) + (FE_B - FE_A)
		         * 
		         * 	[  deltaFE = FE_B - FE_A  ]
		         *  
		         *  
		         *  A: target system
		         *  B: reference system
		         */
		        
		        double sHarmonic = wHarmonic - deltaFE_DSHarmonic;
		        double sTarget =   wTarget + deltaFE_DSHarmonic;
		        
		        DataGroup allYourBase0 = (DataGroup)sim.accumulators[0].getData(sim.dsvo.minDiffLocation());
		        double ratioHarmonicAverage = ((DataDoubleArray)allYourBase0.getData(AccumulatorAverage.StatType.AVERAGE.index)).getData()[1];
		        double deltaFEHarmonic = -Math.log(ratioHarmonicAverage);
		        
		        double biasF = -temp*Math.log(0.5 * SpecialFunctions.erfc( Math.sqrt(sHarmonic)* 
		        		( 1- Math.sqrt( (0.5/sTarget)*Math.log ((0.5/Math.PI)*(idStep-1)*(idStep-1))))));
		        
		        try {
		        	
		        	fileWriterHarm.write(idStep + " " + temp*deltaFE_DSHarmonic + " " + temp*deltaFEHarmonic + " "+ biasF +  "\n");
		        	
		        } catch (IOException e){
		        	
		        }
        	}

        	
        };
        
        IAction outputActionTarget = new IAction(){
        	public void actionPerformed(){

        		long idStep = sim.integrators[1].getStepCount(); // /sim.integrators[1].getActionInterval(sim.accumulatorPumps[1]);
     
		        /*
		         * Histogram
		         */
		        //Target ---> Harmonic
				DataLogger dataLogger1 = new DataLogger();
				DataTableWriter dataTableWriter1 = new DataTableWriter();
				dataLogger1.setFileName(outFileName + "_histDS_TargHarm");
				dataLogger1.setDataSink(dataTableWriter1);
				dataTableWriter1.setIncludeHeader(false);
				dataLogger1.putDataInfo(histogramTargetHarmonic.getDataInfo());
				
				dataLogger1.setWriteInterval(1);
				dataLogger1.setAppending(false); //overwrite data
				dataLogger1.putData(histogramTargetHarmonic.getData());
				
				dataLogger1.closeFile();

		        double[][] omega2 = sim.normalModes.getOmegaSquared(sim.boxTarget);
		        double[] coeffs = sim.normalModes.getWaveVectorFactory().getCoefficients();
		        double AHarmonic = 0;
		        for(int i=0; i<omega2.length; i++) {
		            for(int j=0; j<omega2[0].length; j++) {
		                if (!Double.isInfinite(omega2[i][j])) {
		                    AHarmonic += coeffs[i]*Math.log(omega2[i][j]*coeffs[i]/(temp*Math.PI));
		                }
		            }           
		            
		        }
		
		        int totalCells = 1;
		        for (int i=0; i<dimension; i++) {
		            totalCells *= sim.nCells[i];
		        }
		        int basisSize = sim.basis.getScaledCoordinates().length;
		        double fac = 1;
		        if (totalCells % 2 == 0) {
		            fac = Math.pow(2,dimension);
		        }
		        AHarmonic -= Math.log(Math.pow(2.0, basisSize*dimension*(totalCells - fac)/2.0) / Math.pow(totalCells,0.5*dimension));
		     
		        /*
		         * ratio = Q_target/Q_reference
		         * delta_FE = FE_target- FE_reference = - ln (ratio)
		         * 
		         * FE_target = FE_reference + deltaFE
		         * 
		         * free energy difference = temp * (FE_reference + FE_target)
		         * 
		         *  note: Aharmonic has not been timed temperature yet
		         *        that is why FE_reference * temp
		         */
		        double ratio = sim.dsvo.getDataAsScalar();
		        double deltaFE = -Math.log(ratio);
		       
		        double boltzmannTargetAverage = dataAverageBoltzmannTarget.getData().getValue(AccumulatorAverage.StatType.AVERAGE.index);
		        double deltaFE_DSTarget = -Math.log(boltzmannTargetAverage);
        
		        double wHarmonic = dataAverageHarmonic.getData().getValue(AccumulatorAverage.StatType.AVERAGE.index);
		        double wTarget = dataAverageTarget.getData().getValue(AccumulatorAverage.StatType.AVERAGE.index);
		
		        /*
		         * s_A = (u_B - u_A) - (FE_B - FE_A) 
		         * s_B = (u_A - u_B) + (FE_B - FE_A)
		         * 
		         * 	[  deltaFE = FE_B - FE_A  ]
		         *  
		         *  
		         *  A: target system
		         *  B: reference system
		         */
		        
		        double sHarmonic = wHarmonic - deltaFE_DSTarget;
		        double sTarget =   wTarget + deltaFE_DSTarget;
		        
		        DataGroup allYourBase1 = (DataGroup)sim.accumulators[1].getData(sim.dsvo.minDiffLocation());
		        double ratioTargetAverage = ((DataDoubleArray)allYourBase1.getData(AccumulatorAverage.StatType.AVERAGE.index)).getData()[1];
		        double deltaFETarget = -Math.log(ratioTargetAverage);
		        
		        
		        double biasR = -temp*Math.log(0.5 * SpecialFunctions.erfc( Math.sqrt(sTarget)* 
		        		( 1- Math.sqrt( (0.5/sHarmonic)*Math.log ((0.5/Math.PI)*(idStep-1)*(idStep-1))))));
		        
		        try {
		        	
		        	fileWriterTarg.write(idStep + " " + temp*deltaFE_DSTarget + " " + temp*deltaFETarget+ " "+ biasR +"\n");
		        	
		        } catch (IOException e){
		        	
		        }
        	}

        	
        };
        
        sim.integratorOverlap.addIntervalAction(outputActionOverlap);
        sim.integratorOverlap.setActionInterval(outputActionOverlap, 100); //100
        
        sim.integrators[0].addIntervalAction(outputActionHarmonic);
        sim.integrators[0].setActionInterval(outputActionHarmonic, 10000); //2000
        
        sim.integrators[1].addIntervalAction(outputActionTarget);
        sim.integrators[1].setActionInterval(outputActionTarget, 100000); //8000
        
        sim.activityIntegrate.setMaxSteps(numSteps);
        sim.getController().actionPerformed();
        
        System.out.println("final reference optimal step frequency "+sim.integratorOverlap.getStepFreq0()+" (actual: "+sim.integratorOverlap.getActualStepFreq0()+")");
        
        double[][] omega2 = sim.normalModes.getOmegaSquared(sim.boxTarget);
        double[] coeffs = sim.normalModes.getWaveVectorFactory().getCoefficients();
        double AHarmonic = 0;
        for(int i=0; i<omega2.length; i++) {
            for(int j=0; j<omega2[0].length; j++) {
                if (!Double.isInfinite(omega2[i][j])) {
                    AHarmonic += coeffs[i]*Math.log(omega2[i][j]*coeffs[i]/(temperature*Math.PI));
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
        AHarmonic -= Math.log(Math.pow(2.0, basisSize*D*(totalCells - fac)/2.0) / Math.pow(totalCells,0.5*D));
        System.out.println("Harmonic-reference free energy: "+AHarmonic*temperature);

        double ratio = sim.dsvo.getDataAsScalar();
        double error = sim.dsvo.getError();
        System.out.println("ratio average: "+ratio+", error: "+error);
        System.out.println("free energy difference: "+(-temperature*Math.log(ratio))+" ,error: "+temperature*(error/ratio));
        System.out.println("target free energy: "+temperature*(AHarmonic-Math.log(ratio))+" ,error: "+temperature*(error/ratio));
        System.out.println("target free energy per particle: "+temperature*(AHarmonic-Math.log(ratio))/numMolecules);
        DataGroup allYourBase = (DataGroup)sim.accumulators[0].getData(sim.dsvo.minDiffLocation());
        System.out.println("harmonic ratio average: "+((DataDoubleArray)allYourBase.getData(AccumulatorAverage.StatType.AVERAGE.index)).getData()[1]
                          +" stdev: "+((DataDoubleArray)allYourBase.getData(AccumulatorAverage.StatType.STANDARD_DEVIATION.index)).getData()[1]
                          +" error: "+((DataDoubleArray)allYourBase.getData(AccumulatorAverage.StatType.ERROR.index)).getData()[1]);
        
        allYourBase = (DataGroup)sim.accumulators[1].getData(sim.dsvo.minDiffLocation());
        System.out.println("target ratio average: "+((DataDoubleArray)allYourBase.getData(AccumulatorAverage.StatType.AVERAGE.index)).getData()[1]
                          +" stdev: "+((DataDoubleArray)allYourBase.getData(AccumulatorAverage.StatType.STANDARD_DEVIATION.index)).getData()[1]
                          +" error: "+((DataDoubleArray)allYourBase.getData(AccumulatorAverage.StatType.ERROR.index)).getData()[1]);
        
        System.out.println("Average-Pressure: "+ ((DataGroup)pressureTargetAverage.getData()).getValue(AccumulatorAverage.StatType.AVERAGE.index) +
				 " ,Error-Pressure: "+ ((DataGroup)pressureTargetAverage.getData()).getValue(AccumulatorAverage.StatType.ERROR.index));

        try{
	        fileWriterDSFE.close();
	        fileWriterBen.close();
	        fileWriterHarm.close();
	        fileWriterTarg.close();
        } catch (IOException e){
        	
        }
        
        
    }

    private static final long serialVersionUID = 1L;
    public IntegratorOverlap integratorOverlap;
    public DataSourceVirialOverlap dsvo;
    public IntegratorBox[] integrators;
    public ActivityIntegrate activityIntegrate;
    public IBox boxTarget, boxHarmonic;
    public Boundary boundaryTarget, boundaryHarmonic;
    public int[] nCells;
    public Basis basis;
    public NormalModes normalModes;
    public Primitive primitive;
    public double refPref;
    public AccumulatorVirialOverlapSingleAverage[] accumulators;
    public DataPump[] accumulatorPumps;
    public IEtomicaDataSource[] meters;
    public AccumulatorAverageFixed accumulatorWorkAverage;
    public MCMoveHarmonic move;
    public IntegratorMC integratorTarget;
    public MeterHarmonicEnergy meterHarmonicEnergy;
    public PotentialMasterMonatomic potentialMasterTarget;
    public double latticeEnergy;

    /**
     * Inner class for parameters understood by the HSMD3D constructor
     */
    public static class SimOverlapParam extends ParameterBase {
        public int numMolecules =32;
        public double density = 1256;
        public int exponentN = 12;
        public int D = 3;
        public long numSteps = 1000000;
        public double harmonicFudge = 1;
        public String filename = "DB_FCC_n12_T16";
        public double temperature = 1.6;
    }
}
