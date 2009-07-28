package etomica.normalmode;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;

import etomica.action.IAction;
import etomica.action.activity.ActivityIntegrate;
import etomica.api.IAtomType;
import etomica.api.IBox;
import etomica.box.Box;
import etomica.data.AccumulatorAverage;
import etomica.data.AccumulatorAverageFixed;
import etomica.data.AccumulatorHistogram;
import etomica.data.DataFork;
import etomica.data.DataLogger;
import etomica.data.DataPump;
import etomica.data.DataTableWriter;
import etomica.data.IEtomicaDataSource;
import etomica.data.meter.MeterPotentialEnergy;
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
import etomica.listener.IntegratorListenerAction;
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
 * Simulation to run sampling with the hard sphere potential, but measuring
 * the harmonic potential based on normal mode data from a previous simulation.
 * 
 * @author Andrew Schultz & Tai Tan
 */
public class SimPhaseSpaceOverlapSoftSphere extends Simulation {

    public SimPhaseSpaceOverlapSoftSphere(Space _space, int numAtoms, double density, double temperature, String filename, double harmonicFudge, int exponent) {
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
        double truncationRadius = boundaryTarget.getBoxSize().getX(0) * 0.495;
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
        integratorOverlap = new IntegratorOverlap(new IntegratorBox[]{integratorHarmonic, integratorTarget});
        meterHarmonicEnergy = new MeterHarmonicEnergy(coordinateDefinitionTarget, normalModes);
        
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
        newAccumulator.setBlockSize(100);
        if (accumulatorPumps[iBox] == null) {
            accumulatorPumps[iBox] = new DataPump(meters[iBox],newAccumulator);
            IntegratorListenerAction pumpListener = new IntegratorListenerAction(accumulatorPumps[iBox]);
            integrators[iBox].getEventManager().addListener(pumpListener);
            if (iBox == 1) {
                pumpListener.setInterval(boxTarget.getMoleculeList().getMoleculeCount());
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
        setAccumulator(new AccumulatorVirialOverlapSingleAverage(1,true),0);
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
                setAccumulator(new AccumulatorVirialOverlapSingleAverage(1,true),0);
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

            setAccumulator(new AccumulatorVirialOverlapSingleAverage(41,true),0);
            setAccumulator(new AccumulatorVirialOverlapSingleAverage(41,false),1);
            setRefPref(1,200);
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
            
            setAccumulator(new AccumulatorVirialOverlapSingleAverage(11,true),0);
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
            setAccumulator(new AccumulatorVirialOverlapSingleAverage(1,true),0);
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
     * @see SimPhaseSpaceOverlapSoftSphere.SimOverlapParam
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
        String refFileName = args.length > 0 ? filename+"_ref" : null;

        System.out.println("Running "+(D==1 ? "1D" : (D==3 ? "FCC" : "2D hexagonal")) +" soft sphere overlap simulation");
        System.out.println(numMolecules+" atoms at density "+density+" and temperature "+temperature);
        System.out.println("exponent N: "+ exponentN +" and harmonic fudge: "+harmonicFudge);
        System.out.println((numSteps/1000)+" total steps of 1000");
        System.out.println("output data to "+filename);

        //instantiate simulation
        final SimPhaseSpaceOverlapSoftSphere sim = new SimPhaseSpaceOverlapSoftSphere(Space.getInstance(D), numMolecules, density, temperature, filename, harmonicFudge, exponentN);
        
        //start simulation
        sim.integratorOverlap.setNumSubSteps(1000);
        numSteps /= 1000;

//        StopWatcher stopWatcher = new StopWatcher("stop", sim, filename);
//        IntervalActionAdapter iaa = new IntervalActionAdapter(stopWatcher);
//        iaa.setActionInterval(100);
//        sim.integratorOverlap.addListener(iaa);

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
        IEtomicaDataSource[] boltzmannDirectSampling = new IEtomicaDataSource[1];
        // Work Determination 
        
        // Target --> Harmonic
        MeterWorkTargetPhaseSpace meterWorkTarget =  new MeterWorkTargetPhaseSpace(sim.integratorTarget, sim.meterHarmonicEnergy);
        meterWorkTarget.setLatticeEnergy(sim.latticeEnergy);
        workMeters[1] = meterWorkTarget;
        
        DataFork dataForkTarget = new DataFork();
        DataPump pumpTarget = new DataPump(workMeters[1], dataForkTarget);
        
        final AccumulatorAverageFixed dataAverageTarget = new AccumulatorAverageFixed();
        
        dataForkTarget.addDataSink(dataAverageTarget);
        IntegratorListenerAction pumpTargetListener = new IntegratorListenerAction(pumpTarget);
        pumpTargetListener.setInterval(1);
        sim.integrators[1].getEventManager().addListener(pumpTargetListener);
        
        //Histogram Target ---> Harmonic
        final AccumulatorHistogram histogramTargetHarmonic = new AccumulatorHistogram(new HistogramSimple(600, new DoubleRange(-150, 450)));
        dataForkTarget.addDataSink(histogramTargetHarmonic);
        
        
     //Direct Sampling - Boltzmann Factor difference average (Harmonic ---> Target)  
        MeterDirectSamplingHarmonic meterDirectSampling = new MeterDirectSamplingHarmonic(sim.move, sim.potentialMasterTarget);
        meterDirectSampling.setTemperature(temperature);
        meterDirectSampling.setLatticeEnergy(sim.latticeEnergy);
        
        boltzmannDirectSampling[0] = meterDirectSampling;
        
        final AccumulatorAverageFixed dataAverageBoltzmann = new AccumulatorAverageFixed();
        DataPump pumpBoltzmann = new DataPump(boltzmannDirectSampling[0], dataAverageBoltzmann);
        
        IntegratorListenerAction pumpBoltzmannListener = new IntegratorListenerAction(pumpBoltzmann);
        pumpBoltzmannListener.setInterval(1);
        sim.integrators[0].getEventManager().addListener(pumpBoltzmannListener);
        
        // Harmonic --> Target
        MeterWorkHarmonicPhaseSpace meterWorkHarmonic = new MeterWorkHarmonicPhaseSpace(sim.move, sim.potentialMasterTarget);
        meterWorkHarmonic.setTemperature(temperature);
        meterWorkHarmonic.setLatticeEnergy(sim.latticeEnergy);
        workMeters[0] = meterWorkHarmonic;
        
        DataFork dataForkHarmonic = new DataFork();
        DataPump pumpHarmonic = new DataPump(workMeters[0], dataForkHarmonic);
        
        final AccumulatorAverageFixed dataAverageHarmonic = new AccumulatorAverageFixed();
        
        dataForkHarmonic.addDataSink(dataAverageHarmonic);
        IntegratorListenerAction pumpHarmonicListener = new IntegratorListenerAction(pumpHarmonic);
        pumpHarmonicListener.setInterval(1);
        sim.integrators[0].getEventManager().addListener(pumpHarmonicListener);
        
        /*
         * Perturbation 
         */
      // Target ---> Bennet's Overlap
        MeterWorkTargetBennet meterWorkTargetBennet = new MeterWorkTargetBennet(sim.integratorTarget, sim.meterHarmonicEnergy,sim.refPref);
        meterWorkTargetBennet.setLatticeEnergy(sim.latticeEnergy);
        workBennets[1] = meterWorkTargetBennet;
        
        DataFork dataForkTargetBennet = new DataFork();
        DataPump pumpTargetBennet = new DataPump(workBennets[1], dataForkTargetBennet);
        
        final AccumulatorAverageFixed dataAverageTargetBennet = new AccumulatorAverageFixed();
        
        dataForkTargetBennet.addDataSink(dataAverageTargetBennet);
        IntegratorListenerAction pumpTargetBennetListener = new IntegratorListenerAction(pumpTargetBennet);
        pumpTargetBennetListener.setInterval(numMolecules*2);
        sim.integrators[1].getEventManager().addListener(pumpTargetBennetListener);
        
        //Histogram Target
        final AccumulatorHistogram histogramTarget = new AccumulatorHistogram(new HistogramSimple(600,new DoubleRange(-150,450)));
        dataForkTargetBennet.addDataSink(histogramTarget);
        
        
      // Harmonic ---> Bennet's Overlap
        MeterWorkHarmonicBennet meterWorkHarmonicBennet = new MeterWorkHarmonicBennet(sim.move, sim.potentialMasterTarget, sim.refPref);
        meterWorkHarmonicBennet.setTemperature(temperature);
        meterWorkHarmonicBennet.setLatticeEnergy(sim.latticeEnergy);
        workBennets[0] = meterWorkHarmonicBennet;
        
        DataFork dataForkHarmonicBennet = new DataFork();
        DataPump pumpHarmonicBennet = new DataPump(workBennets[0], dataForkHarmonicBennet);
        
        final AccumulatorAverageFixed dataAverageHarmonicBennet = new AccumulatorAverageFixed();
        
        dataForkHarmonicBennet.addDataSink(dataAverageHarmonicBennet);
        IntegratorListenerAction pumpHarmonicBennetListener = new IntegratorListenerAction(pumpHarmonicBennet);
        pumpHarmonicBennetListener.setInterval(1);
        sim.integrators[0].getEventManager().addListener(pumpHarmonicBennetListener);
        
        //Histogram Harmonic
        final AccumulatorHistogram histogramHarmonic = new AccumulatorHistogram(new HistogramSimple(600, new DoubleRange(-150,450)));
        dataForkHarmonicBennet.addDataSink(histogramHarmonic);
        
        FileWriter fileWriter, fileWriterBennet;
               
        try{
        	fileWriter = new FileWriter(filename+"_FE");
        	fileWriterBennet = new FileWriter(filename+"_Ben");
        }catch(IOException e){
        	fileWriter = null;
        	fileWriterBennet = null;
        }
        
        final String outFileName = filename;
        final double temp = temperature;
        final int dimension = D;
        final long steps = numSteps;
        final FileWriter fileWriterFE = fileWriter;
        final FileWriter fileWriterBen = fileWriterBennet;
        
        IAction outputAction = new IAction(){
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
		        
				//Harmonic
				DataLogger dataLogger2 = new DataLogger();
				DataTableWriter dataTableWriter2 = new DataTableWriter();
				dataLogger2.setFileName(outFileName + "_hist_HarmBenn");
				dataLogger2.setDataSink(dataTableWriter2);
				dataTableWriter2.setIncludeHeader(false);
				dataLogger2.putDataInfo(histogramHarmonic.getDataInfo());
				
				dataLogger2.setWriteInterval(1);
				dataLogger2.setAppending(false); //overwrite data
				dataLogger2.putData(histogramHarmonic.getData());
				dataLogger2.closeFile();
			
		        
		        //Target
				DataLogger dataLogger3 = new DataLogger();
				DataTableWriter dataTableWriter3 = new DataTableWriter();
				dataLogger3.setFileName(outFileName + "_hist_TargBenn");
				dataLogger3.setDataSink(dataTableWriter3);
				dataTableWriter3.setIncludeHeader(false);
				dataLogger3.putDataInfo(histogramTarget.getDataInfo());
				
				dataLogger3.setWriteInterval(1);
				dataLogger3.setAppending(false); //overwrite data
				dataLogger3.putData(histogramTarget.getData());
				dataLogger3.closeFile();
				
		        /*
		         * 
		         */
		        
		        
		        
		        
		        System.out.println("\n*****************************************************************");
		        System.out.println("**********              "+ idStep*1000 + "            ***************");
		        System.out.println("*****************************************************************");
		        System.out.println("\nfinal reference optimal step frequency "+sim.integratorOverlap.getStepFreq0()+" (actual: "+sim.integratorOverlap.getActualStepFreq0()+")");
		        
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
		        System.out.println("Harmonic-reference free energy: "+AHarmonic*temp);
		        
		        
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
		        
		        double boltzmannAverage = dataAverageBoltzmann.getData().getValue(AccumulatorAverage.StatType.AVERAGE.index);
		        double boltzmannError = dataAverageBoltzmann.getData().getValue(AccumulatorAverage.StatType.ERROR.index);
		        
		        double deltaFE_DirectSampling = -Math.log(boltzmannAverage);
		        System.out.println("boltzmann: " + boltzmannAverage);
		        System.out.println("ratio average: "+ratio+", error: "+error);
		        System.out.println("free energy difference (Bennett): "+ (temp*deltaFE) +" ,error: "+ errorFE);
		        System.out.println("free energy difference (Direct Sampling): "+ (temp*deltaFE_DirectSampling) 
		        											+ " ,error: " + boltzmannError/boltzmannAverage);
		        
		        System.out.println("target free energy: "+ targetFE);
		        
		        

		                
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
		        
		        System.out.println("wHarmonic: "+ wHarmonic + 
		        		           " ,error: "+ eHarmonic);
		        System.out.println("wTarget: " + wTarget +
		        		           " ,error: "+ eTarget);
		        System.out.println("beta*deltaFE (Bennett): " + deltaFE);
		        System.out.println("beta*deltaFE (Direct Sampling): " + deltaFE_DirectSampling);
		        
		        System.out.println("\nHarmonic entropy: " + sHarmonic + " ,error: "+ er_sHarmonic);
		        System.out.println("Target entropy: " + sTarget + " ,error: "+ er_sTarget);
		        System.out.println("ratio entropy (harmonic-to-target) "+ sHarmonic/sTarget);
		         
		        double mNumStepsHarmonic = steps*1000*sim.integratorOverlap.getActualStepFreq0();
		        double mNumStepsTarget = steps*1000*(1-sim.integratorOverlap.getActualStepFreq0());
		        double pi_harmonic = Math.sqrt((sHarmonic/sTarget)*Math.log((0.5/Math.PI)*mNumStepsHarmonic*mNumStepsHarmonic))-Math.sqrt(2*sHarmonic);
		        double pi_target = Math.sqrt((sTarget/sHarmonic)*Math.log((0.5/Math.PI)*mNumStepsTarget*mNumStepsTarget))-Math.sqrt(2*sTarget);
		        
		        System.out.println("PI Harmonic: " + pi_harmonic + " ,Msteps: " + mNumStepsHarmonic);
		        System.out.println("PI Target: " + pi_target + " ,Msteps: " + mNumStepsTarget);
		        System.out.println();
		        
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
		        
		        System.out.println("harmonic ratio average: "+ ratioHarmonicAverage
		                          +" stdev: "+((DataDoubleArray)allYourBase0.getData(AccumulatorAverage.StatType.STANDARD_DEVIATION.index)).getData()[1]
		                          +" error: "+ ratioHarmonicError);
		        
		        DataGroup allYourBase1 = (DataGroup)sim.accumulators[1].getData(sim.dsvo.minDiffLocation());
		        double ratioTargetAverage = ((DataDoubleArray)allYourBase1.getData(AccumulatorAverage.StatType.AVERAGE.index)).getData()[1];
		        double ratioTargetError = ((DataDoubleArray)allYourBase1.getData(AccumulatorAverage.StatType.ERROR.index)).getData()[1];
		        double deltaFETarget = -Math.log(ratioTargetAverage);
		        
		        System.out.println("target ratio average: "+ratioTargetAverage
		                          +" stdev: "+((DataDoubleArray)allYourBase1.getData(AccumulatorAverage.StatType.STANDARD_DEVIATION.index)).getData()[1]
		                          +" error: "+ratioTargetError);
		        
		        double sHarmonicBennet = wHarmonicBennet - deltaFEHarmonic;
		        double sTargetBennet = wTargetBennet - deltaFETarget;
		        
		        double er_sHarmonicBennet = Math.sqrt(eHarmonicBennet*eHarmonicBennet + (ratioHarmonicError/ratioHarmonicAverage)*(ratioHarmonicError/ratioHarmonicAverage));
		        double er_sTargetBennet = Math.sqrt(eTargetBennet*eTargetBennet + (ratioTargetError/ratioTargetAverage)*(ratioTargetError/ratioTargetAverage));
		        
		        System.out.println("\nwHarmonicBennet: "+ wHarmonicBennet + 
				           " ,error: "+ eHarmonicBennet);
		        System.out.println("wTargetBennet: " + wTargetBennet +
				           " ,error: "+ eTargetBennet);
		        System.out.println("deltaFE Harmonic: "+ deltaFEHarmonic);
		        System.out.println("deltaFE Target: "+ deltaFETarget);
		        
		        System.out.println("\nHarmonic entropy (perturbed into Bennet): " + sHarmonicBennet + 
		        		" ,error: "+ er_sHarmonicBennet);
		        System.out.println("Target entropy (perturbed into Bennet): " + sTargetBennet + 
		        		" ,error: "+ er_sTargetBennet);
		        
		        double pi_harmonicBennet = Math.sqrt((sHarmonicBennet/sTargetBennet)*Math.log((0.5/Math.PI)*mNumStepsHarmonic*mNumStepsHarmonic))-Math.sqrt(2*sHarmonicBennet);
		        double pi_targetBennet = Math.sqrt((sTargetBennet/sHarmonicBennet)*Math.log((0.5/Math.PI)*mNumStepsTarget*mNumStepsTarget))-Math.sqrt(2*sTargetBennet);
		        
		        System.out.println("PI Harmonic: " + pi_harmonicBennet + " ,Msteps: " + mNumStepsHarmonic);
		        System.out.println("PI Target: " + pi_targetBennet + " ,Msteps: " + mNumStepsTarget);
		        System.out.println();
		        
		        try {
		        	
		        	fileWriterFE.write(idStep*1000 + " " + temp*deltaFE + " " + temp*deltaFE_DirectSampling  
		        							  + " " + targetFE + " "+ errorFE + " "+ wHarmonic + " " + wTarget
		        			                  + " " + sHarmonic+ " "+ sTarget + " "+ pi_harmonic + " " + pi_target +"\n");
		        	
		        	fileWriterBen.write(idStep*1000 + " " + deltaFEHarmonic + " " + deltaFETarget + " " + wHarmonicBennet + " "+ wTargetBennet
		        			                      + " " + sHarmonicBennet + " " + sTargetBennet + " " + pi_harmonicBennet + " "+ pi_targetBennet + "\n");
		        	
		        } catch (IOException e){
		        	
		        }
        	}
        }; 
        
        IntegratorListenerAction outputActionListener = new IntegratorListenerAction(outputAction);
        outputActionListener.setInterval(10);
        sim.integratorOverlap.getEventManager().addListener(outputActionListener);
        
        sim.activityIntegrate.setMaxSteps(numSteps);
        sim.getController().actionPerformed();
        
        try{
	        fileWriterFE.close();
	        fileWriterBen.close();
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
        public String filename = "CB_FCC_n12_T01";
        public double temperature = 0.1;
    }
}
