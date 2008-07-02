package etomica.normalmode;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;

import etomica.action.activity.ActivityIntegrate;
import etomica.api.IAtomTypeLeaf;
import etomica.api.IBox;
import etomica.box.Box;
import etomica.data.AccumulatorAverage;
import etomica.data.AccumulatorAverageFixed;
import etomica.data.AccumulatorHistogram;
import etomica.data.DataLogger;
import etomica.data.DataPump;
import etomica.data.DataSource;
import etomica.data.DataTableWriter;
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
import etomica.potential.P2SoftSphere;
import etomica.potential.P2SoftSphericalTruncated;
import etomica.potential.Potential2SoftSpherical;
import etomica.potential.PotentialMasterMonatomic;
import etomica.simulation.Simulation;
import etomica.space.Boundary;
import etomica.space.BoundaryRectangularPeriodic;
import etomica.space.Space;
import etomica.species.SpeciesSpheresMono;
import etomica.util.DoubleRange;
import etomica.util.HistogramExpanding;
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

        potentialMasterTarget = new PotentialMasterMonatomic(this, space);
        integrators = new IntegratorBox[2];
        accumulatorPumps = new DataPump[2];
        meters = new DataSource[2];
        accumulators = new AccumulatorVirialOverlapSingleAverage[2];

        SpeciesSpheresMono species = new SpeciesSpheresMono(this, space);
        getSpeciesManager().addSpecies(species);

        // TARGET
        boxTarget = new Box(this, space);
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
            boundaryTarget = new BoundaryRectangularPeriodic(space, getRandom(), numAtoms/density);
            nCells = new int[]{numAtoms};
            basis = new BasisMonatomic(space);
        } else {
            double L = Math.pow(4.0/density, 1.0/3.0);
            primitive = new PrimitiveCubic(space, L);
            int n = (int)Math.round(Math.pow(numAtoms/4, 1.0/3.0));
            nCells = new int[]{n,n,n};
            boundaryTarget = new BoundaryRectangularPeriodic(space, random, n * L);
            basis = new BasisCubicFcc();
        }
        boxTarget.setBoundary(boundaryTarget);

        CoordinateDefinitionLeaf coordinateDefinitionTarget = new CoordinateDefinitionLeaf(this, boxTarget, primitive, basis, space);
        coordinateDefinitionTarget.initializeCoordinates(nCells);

        Potential2SoftSpherical potential = new P2SoftSphere(space, 1.0, 1.0, exponent);
        double truncationRadius = boundaryTarget.getDimensions().x(0) * 0.495;
        P2SoftSphericalTruncated pTruncated = new P2SoftSphericalTruncated(space, potential, truncationRadius);
        IAtomTypeLeaf sphereType = species.getLeafType();
        potentialMasterTarget.addPotential(pTruncated, new IAtomTypeLeaf[] { sphereType, sphereType });
        atomMove.setPotential(pTruncated);

        integratorTarget.setBox(boxTarget);

        potentialMasterTarget.lrcMaster().setEnabled(false);
        MeterPotentialEnergy meterPE = new MeterPotentialEnergy(potentialMasterTarget);
        meterPE.setBox(boxTarget);
        latticeEnergy = meterPE.getDataAsScalar();
        
    
        // HARMONIC
        boundaryHarmonic = new BoundaryRectangularPeriodic(random, space);
        boxHarmonic = new Box(boundaryHarmonic, space);
        addBox(boxHarmonic);
        boxHarmonic.setNMolecules(species, numAtoms);

        IntegratorMC integratorHarmonic = new IntegratorMC(potentialMasterTarget, random, 1.0);

        move = new MCMoveHarmonic(getRandom());
        integratorHarmonic.getMoveManager().addMCMove(move);
        integrators[0] = integratorHarmonic;
        
        if (space.D() == 1) {
            boundaryHarmonic = new BoundaryRectangularPeriodic(space, getRandom(), numAtoms/density);
        } else {
            double L = Math.pow(4.0/density, 1.0/3.0);
            int n = (int)Math.round(Math.pow(numAtoms/4, 1.0/3.0));
            boundaryHarmonic = new BoundaryRectangularPeriodic(space, random, n * L);
        }
        boxHarmonic.setBoundary(boundaryHarmonic);

        CoordinateDefinitionLeaf coordinateDefinitionHarmonic = new CoordinateDefinitionLeaf(this, boxHarmonic, primitive, basis, space);
        coordinateDefinitionHarmonic.initializeCoordinates(nCells);
        
        normalModes = new NormalModesFromFile(filename, space.D());
        normalModes.setHarmonicFudge(harmonicFudge);
        
        /*
         * nuke this line if it is overlap between DB and harmonic
         */
        //normalModes.setTemperature(temperature);
        
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
        newAccumulator.setBlockSize(100);
        if (accumulatorPumps[iBox] == null) {
            accumulatorPumps[iBox] = new DataPump(meters[iBox],newAccumulator);
            integrators[iBox].addIntervalAction(accumulatorPumps[iBox]);
            if (iBox == 1) {
                integrators[iBox].setActionInterval(accumulatorPumps[iBox], boxTarget.getMoleculeList().getAtomCount());
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
        SimPhaseSpaceOverlapSoftSphere sim = new SimPhaseSpaceOverlapSoftSphere(Space.getInstance(D), numMolecules, density, temperature, filename, harmonicFudge, exponentN);
        
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

//

        DataSource[] workMeters = new DataSource[2];
        DataSource[] workBennets = new DataSource[2];
        
        // Work Determination 
        
        // Target --> Harmonic
        MeterWorkTargetPhaseSpace meterWorkTarget =  new MeterWorkTargetPhaseSpace(sim.integratorTarget, sim.meterHarmonicEnergy);
        meterWorkTarget.setLatticeEnergy(sim.latticeEnergy);
        workMeters[1] = meterWorkTarget;
        
        AccumulatorAverageFixed dataAverageTarget = new AccumulatorAverageFixed();
        DataPump pumpTarget = new DataPump(workMeters[1], dataAverageTarget);
        sim.integrators[1].addIntervalAction(pumpTarget);
        sim.integrators[1].setActionInterval(pumpTarget, 1);
        
        // Harmonic --> Target
        MeterWorkHarmonicPhaseSpace meterWorkHarmonic = new MeterWorkHarmonicPhaseSpace(sim.move, sim.potentialMasterTarget);
        meterWorkHarmonic.setTemperature(temperature);
        meterWorkHarmonic.setLatticeEnergy(sim.latticeEnergy);
        workMeters[0] = meterWorkHarmonic;
        
        AccumulatorAverageFixed dataAverageHarmonic = new AccumulatorAverageFixed();
        DataPump pumpHarmonic = new DataPump(workMeters[0], dataAverageHarmonic);
        sim.integrators[0].addIntervalAction(pumpHarmonic);
        sim.integrators[0].setActionInterval(pumpHarmonic, 1);
       
        /*
         * Perturbation 
         */
      // Target ---> Bennet's Overlap
        MeterWorkTargetBennet meterWorkTargetBennet = new MeterWorkTargetBennet(sim.integratorTarget, sim.meterHarmonicEnergy,sim.refPref);
        meterWorkTargetBennet.setLatticeEnergy(sim.latticeEnergy);
        workBennets[1] = meterWorkTargetBennet;
        
        AccumulatorAverageFixed dataAverageTargetBennet = new AccumulatorAverageFixed();
        DataPump pumpTargetBennet = new DataPump(workBennets[1], dataAverageTargetBennet);
        sim.integrators[1].addIntervalAction(pumpTargetBennet);
        sim.integrators[1].setActionInterval(pumpTargetBennet, numMolecules*2);
        
        //Histogram Target
        AccumulatorHistogram histogramTarget = new AccumulatorHistogram(new HistogramSimple(600,new DoubleRange(-150,450)));
        DataPump pumpHistogramTarget = new DataPump(workBennets[1], histogramTarget);
        sim.integrators[1].addIntervalAction(pumpHistogramTarget);
        sim.integrators[1].setActionInterval(pumpHistogramTarget, numMolecules*2); //the interval is based on the system size
        
      // Harmonic ---> Bennet's Overlap
        MeterWorkHarmonicBennet meterWorkHarmonicBennet = new MeterWorkHarmonicBennet(sim.move, sim.potentialMasterTarget, sim.refPref);
        meterWorkHarmonicBennet.setTemperature(temperature);
        meterWorkHarmonicBennet.setLatticeEnergy(sim.latticeEnergy);
        workBennets[0] = meterWorkHarmonicBennet;
        
        AccumulatorAverageFixed dataAverageHarmonicBennet = new AccumulatorAverageFixed();
        DataPump pumpHarmonicBennet = new DataPump(workBennets[0], dataAverageHarmonicBennet);
        sim.integrators[0].addIntervalAction(pumpHarmonicBennet);
        sim.integrators[0].setActionInterval(pumpHarmonicBennet, 1);
        
        //Histogram Harmonic
        AccumulatorHistogram histogramHarmonic = new AccumulatorHistogram(new HistogramSimple(600, new DoubleRange(-150,450)));
        DataPump pumpHistogramHarmonic = new DataPump(workBennets[0], histogramHarmonic);
        sim.integrators[0].addIntervalAction(pumpHistogramHarmonic);
        sim.integrators[0].setActionInterval(pumpHistogramHarmonic, 1);
        
        sim.activityIntegrate.setMaxSteps(numSteps);
        sim.getController().actionPerformed();
        
        /*
         * Histogram
         */
        //Target
		DataLogger dataLogger = new DataLogger();
		DataTableWriter dataTableWriter = new DataTableWriter();
		dataLogger.setFileName(filename + "_hist_TargBenn");
		dataLogger.setDataSink(dataTableWriter);
		dataTableWriter.setIncludeHeader(false);
		dataLogger.putDataInfo(histogramTarget.getDataInfo());
		
		dataLogger.setWriteInterval(1);
		dataLogger.putData(histogramTarget.getData());
		dataLogger.closeFile();
        
		//Harmonic
        dataLogger.setFileName(filename + "_hist_HarmBenn");
        dataLogger.putData(histogramHarmonic.getData());
        dataLogger.closeFile();
        /*
         * 
         */
        
        

        System.out.println("\nfinal reference optimal step frequency "+sim.integratorOverlap.getStepFreq0()+" (actual: "+sim.integratorOverlap.getActualStepFreq0()+")");
        
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
        
        
        /*
         * ratio = Q_ref/Q_target
         * delta_FE = FE_reference - FE_target = ln (ratio)
         * 
         * FE_target = FE_reference - ln(Q_ref/Q_target)
         * 
         */
        double ratio = sim.dsvo.getDataAsScalar();
        double error = sim.dsvo.getError();
        double deltaFE = Math.log(ratio);
        
        System.out.println("ratio average: "+ratio+", error: "+error);
        System.out.println("free energy difference: "+ (-temperature*deltaFE) +", error: "+temperature*(error/ratio));
        System.out.println("target free energy: "+temperature*(AHarmonic-deltaFE));
                
        double wHarmonic = dataAverageHarmonic.getData().getValue(AccumulatorAverage.StatType.AVERAGE.index);
        double wTarget = dataAverageTarget.getData().getValue(AccumulatorAverage.StatType.AVERAGE.index);

        double eHarmonic = dataAverageHarmonic.getData().getValue(AccumulatorAverage.StatType.ERROR.index);
        double eTarget = dataAverageTarget.getData().getValue(AccumulatorAverage.StatType.ERROR.index);
        
        /*
         * s_A = (u_B - u_A) - (FE_B - FE_A) 
         * s_B = (u_A - u_B) - (FE_A - FE_B)
         * 
         * 	[  deltaFE = FE_A - FE_B  ]
         *  
         *  
         *  A: reference system
         *  B: target system
         */
        
        double sHarmonic = wHarmonic + deltaFE;
        double sTarget =   wTarget - deltaFE;
        
        double er_sHarmonic = Math.sqrt(eHarmonic*eHarmonic + (error/ratio)*(error/ratio));
        double er_sTarget = Math.sqrt(eTarget*eTarget + (error/ratio)*(error/ratio));
        
        System.out.println("wHarmonic: "+ wHarmonic + 
        		           ", error: "+ eHarmonic);
        System.out.println("wTarget: " + wTarget +
        		           ", error: "+ eTarget);
        System.out.println("deltaFE: " + deltaFE);
        
        System.out.println("\nHarmonic entropy: " + sHarmonic + ", error: "+ er_sHarmonic);
        System.out.println("Target entropy: " + sTarget + ", error: "+ er_sTarget);
        System.out.println("ratio entropy (harmonic-to-target) "+ sHarmonic/sTarget);
         
        double mNumStepsHarmonic = numSteps*1000*sim.integratorOverlap.getActualStepFreq0();
        double mNumStepsTarget = numSteps*1000*(1-sim.integratorOverlap.getActualStepFreq0());
        double pi_harmonic = Math.sqrt((sHarmonic/sTarget)*Math.log((0.5/Math.PI)*mNumStepsHarmonic*mNumStepsHarmonic))-Math.sqrt(2*sHarmonic);
        double pi_target = Math.sqrt((sTarget/sHarmonic)*Math.log((0.5/Math.PI)*mNumStepsTarget*mNumStepsTarget))-Math.sqrt(2*sTarget);
        
        System.out.println("PI Harmonic: " + pi_harmonic + ", Msteps: " + mNumStepsHarmonic);
        System.out.println("PI Target: " + pi_target + ", Msteps: " + mNumStepsTarget);
        System.out.println();
        
        /*
         *     Q0/Q1 = ( Q0/Q_overlap ) / ( Q1/Q_overlap )    
         *         
         *         
         *    s_A = (u_B - u_A) + (FE_A - FE_B) 
         * 
         *   	[  deltaFE_harmonic = FE_A - FE_B  ]
         *      [  deltaFE_target   = FE_A - FE_B  ] 
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
        double deltaFEHarmonic = Math.log(ratioHarmonicAverage);
        
        System.out.println("harmonic ratio average: "+ ratioHarmonicAverage
                          +" stdev: "+((DataDoubleArray)allYourBase0.getData(AccumulatorAverage.StatType.STANDARD_DEVIATION.index)).getData()[1]
                          +" error: "+ ratioHarmonicError);
        
        DataGroup allYourBase1 = (DataGroup)sim.accumulators[1].getData(sim.dsvo.minDiffLocation());
        double ratioTargetAverage = ((DataDoubleArray)allYourBase1.getData(AccumulatorAverage.StatType.AVERAGE.index)).getData()[1];
        double ratioTargetError = ((DataDoubleArray)allYourBase1.getData(AccumulatorAverage.StatType.ERROR.index)).getData()[1];
        double deltaFETarget = Math.log(ratioTargetAverage);
        
        System.out.println("target ratio average: "+ratioTargetAverage
                          +" stdev: "+((DataDoubleArray)allYourBase1.getData(AccumulatorAverage.StatType.STANDARD_DEVIATION.index)).getData()[1]
                          +" error: "+ratioTargetError);
        
        double sHarmonicBennet = wHarmonicBennet + deltaFEHarmonic;
        double sTargetBennet = wTargetBennet + deltaFETarget;
        
        double er_sHarmonicBennet = Math.sqrt(eHarmonicBennet*eHarmonicBennet + (ratioHarmonicError/ratioHarmonicAverage)*(ratioHarmonicError/ratioHarmonicAverage));
        double er_sTargetBennet = Math.sqrt(eTargetBennet*eTargetBennet + (ratioTargetError/ratioTargetAverage)*(ratioTargetError/ratioTargetAverage));
        
        System.out.println("\nwHarmonicBennet: "+ wHarmonicBennet + 
		           ", error: "+ eHarmonicBennet);
        System.out.println("wTargetBennet: " + wTargetBennet +
		           ", error: "+ eTargetBennet);
        System.out.println("deltaFE Harmonic: "+ deltaFEHarmonic);
        System.out.println("deltaFE Target: "+ deltaFETarget);
        
        System.out.println("\nHarmonic entropy (perturbed into Bennet): " + sHarmonicBennet + 
        		", error: "+ er_sHarmonicBennet);
        System.out.println("Target entropy (perturbed into Bennet): " + sTargetBennet + 
        		", error: "+ er_sTargetBennet);
        
        double pi_harmonicBennet = Math.sqrt((sHarmonicBennet/sTargetBennet)*Math.log((0.5/Math.PI)*mNumStepsHarmonic*mNumStepsHarmonic))-Math.sqrt(2*sHarmonicBennet);
        double pi_targetBennet = Math.sqrt((sTargetBennet/sHarmonicBennet)*Math.log((0.5/Math.PI)*mNumStepsTarget*mNumStepsTarget))-Math.sqrt(2*sTargetBennet);
        
        System.out.println("PI Harmonic: " + pi_harmonicBennet + ", Msteps: " + mNumStepsHarmonic);
        System.out.println("PI Target: " + pi_targetBennet + ", Msteps: " + mNumStepsTarget);
        System.out.println();
        
        
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
    public DataSource[] meters;
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
        public long numSteps = 5000000;
        public double harmonicFudge = 1;
        public String filename = "DB_FCC_n12_T10";
        public double temperature = 1.0;
    }
}
