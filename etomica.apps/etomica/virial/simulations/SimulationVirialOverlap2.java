package etomica.virial.simulations;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;

import etomica.action.activity.ActivityIntegrate;
import etomica.api.ISpecies;
import etomica.data.DataPumpListener;
import etomica.integrator.IntegratorMC;
import etomica.integrator.mcmove.MCMoveBoxStep;
import etomica.integrator.mcmove.MCMoveManager;
import etomica.overlap.IntegratorOverlap;
import etomica.potential.PotentialMaster;
import etomica.simulation.Simulation;
import etomica.space.ISpace;
import etomica.species.SpeciesSpheresMono;
import etomica.species.SpeciesSpheresRotating;
import etomica.virial.BoxCluster;
import etomica.virial.ClusterAbstract;
import etomica.virial.ClusterWeight;
import etomica.virial.ClusterWeightAbs;
import etomica.virial.ConfigurationCluster;
import etomica.virial.MCMoveClusterAngleBend;
import etomica.virial.MCMoveClusterAtomMulti;
import etomica.virial.MCMoveClusterAtomRotateMulti;
import etomica.virial.MCMoveClusterMoleculeMulti;
import etomica.virial.MCMoveClusterRotateMoleculeMulti;
import etomica.virial.MCMoveClusterWiggleMulti;
import etomica.virial.MeterVirial;
import etomica.virial.overlap.AccumulatorVirialOverlapSingleAverage;
import etomica.virial.overlap.DataSourceVirialOverlap;

/**
 * Simulation implementing the overlap-sampling approach to evaluating a cluster
 * diagram.
 * 
 * @author kofke, Andrew Schultz
 */
public class SimulationVirialOverlap2 extends Simulation {

	/*
	 * If this constructor is used to instantiate the simulation, then doWiggle is set to false, and 
	 * ClusterAbstract[] is set to {refCluster,targetCluster}
	 */
    public SimulationVirialOverlap2(ISpace aSpace, ISpecies species, 
            double temperature, ClusterAbstract refCluster, ClusterAbstract targetCluster) {
        this(aSpace,species,temperature,refCluster, targetCluster, false);
    }

    public SimulationVirialOverlap2(ISpace aSpace, ISpecies species, 
            double temperature, ClusterAbstract refCluster, ClusterAbstract targetCluster, boolean doWiggle) {
        this(aSpace,species,temperature,new ClusterAbstract[]{refCluster,targetCluster},
                new ClusterWeight[]{ClusterWeightAbs.makeWeightCluster(refCluster),ClusterWeightAbs.makeWeightCluster(targetCluster)}, doWiggle);
    }
    
    public SimulationVirialOverlap2(ISpace aSpace, ISpecies species, 
            double temperature, final ClusterAbstract[] aValueClusters, final ClusterWeight[] aSampleClusters, boolean doWiggle) {
        this(aSpace, new ISpecies[]{species}, new int[]{aValueClusters[0].pointCount()}, temperature, aValueClusters, aSampleClusters, doWiggle);
    }

    public SimulationVirialOverlap2(ISpace aSpace, ISpecies[] species, int[] nMolecules,
            double temperature, final ClusterAbstract[] aValueClusters, final ClusterWeight[] aSampleClusters, boolean doWiggle) {

	    super(aSpace);
		PotentialMaster potentialMaster = new PotentialMaster();
        sampleClusters = aSampleClusters;
        boolean doRotate = false;
        boolean multiAtomic = false;
        for (int i=0; i<species.length; i++) {
            addSpecies(species[i]);
            if (!(species[i] instanceof SpeciesSpheresMono)) {
                doRotate = true;
            }
            if (!(species[i] instanceof SpeciesSpheresMono || species[i] instanceof SpeciesSpheresRotating)) {
                multiAtomic = true;
            }
        }
        accumulators = new AccumulatorVirialOverlapSingleAverage[sampleClusters.length];
        accumulatorPumps = new DataPumpListener[sampleClusters.length];
        box = new BoxCluster[sampleClusters.length];
        integrators = new IntegratorMC[sampleClusters.length];
        meters = new MeterVirial[sampleClusters.length];
        mcMoveTranslate = new MCMoveBoxStep[sampleClusters.length];
        if (doRotate) {
            mcMoveRotate = new MCMoveBoxStep[sampleClusters.length];
        }
        if (doWiggle) {
            mcMoveWiggle = new MCMoveBoxStep[sampleClusters.length];
        }
        
        blockSize = 1000;
        
        for (int iBox=0; iBox<sampleClusters.length; iBox++) {
            // integrator for iBox samples based on iBox cluster
            box[iBox] = new BoxCluster(sampleClusters[iBox], space);
            addBox(box[iBox]);
            for (int i=0; i<species.length; i++) {
                box[iBox].setNMolecules(species[i], nMolecules[i]);
            }
            
            integrators[iBox] = new IntegratorMC(this, potentialMaster);
            integrators[iBox].setTemperature(temperature);
            integrators[iBox].setBox(box[iBox]);
            integrators[iBox].getMoveManager().setEquilibrating(true);
            
            MCMoveManager moveManager = integrators[iBox].getMoveManager();
            
            if (!multiAtomic) {
                mcMoveTranslate[iBox] = new MCMoveClusterAtomMulti(random, space);
                moveManager.addMCMove(mcMoveTranslate[iBox]);
                
                if (doRotate) {
                    mcMoveRotate[iBox] = new MCMoveClusterAtomRotateMulti(random, space, aValueClusters[0].pointCount()-1);
                    moveManager.addMCMove(mcMoveRotate[iBox]);
                }
            }
            else {
                mcMoveRotate[iBox] = new MCMoveClusterRotateMoleculeMulti(getRandom(), space);
                mcMoveRotate[iBox].setStepSize(Math.PI);
                moveManager.addMCMove(mcMoveRotate[iBox]);
                mcMoveTranslate[iBox] = new MCMoveClusterMoleculeMulti(this, space);
                moveManager.addMCMove(mcMoveTranslate[iBox]);
                if (doWiggle) {
                    // we can use the bending move if none of the molecules has more than 3 atoms
                    boolean doBend = true;
                    for (int i=0; i<species.length; i++) {
                        if (box[iBox].getMoleculeList(species[i]).getMolecule(0).getChildList().getAtomCount() > 3) {
                            doBend = false;
                        }
                    }
                    if (doBend) {
                        mcMoveWiggle[iBox] = new MCMoveClusterAngleBend(potentialMaster, random, 0.5, space);
                    }
                    else {
                        mcMoveWiggle[iBox] = new MCMoveClusterWiggleMulti(this, potentialMaster, aValueClusters[0].pointCount(), space);
                    }
                    moveManager.addMCMove(mcMoveWiggle[iBox]);
                }
            }
            
            ConfigurationCluster configuration = new ConfigurationCluster(space);
            configuration.initializeCoordinates(box[iBox]);
            MeterVirial meter = new MeterVirial(new ClusterAbstract[]{aValueClusters[iBox],aSampleClusters[1-iBox].makeCopy()});
            setMeter(meter,iBox);
            AccumulatorVirialOverlapSingleAverage acc = new AccumulatorVirialOverlapSingleAverage(11, iBox==0);
            setAccumulator(acc,iBox);
              
        }
        
        setRefPref(1,5);
        integratorOS = new IntegratorOverlap(integrators);
        integratorOS.setNumSubSteps(1000);
        integratorOS.setEventInterval(1);
        ai = new ActivityIntegrate(integratorOS);
        getController().addAction(ai);
		
		dsvo = new DataSourceVirialOverlap(accumulators[0],accumulators[1]);
        integratorOS.setReferenceFracSource(dsvo);
	}

    public void setRefPref(double refPrefCenter, double span) {
        refPref = refPrefCenter;
        accumulators[0].setBennetParam(refPrefCenter,span);
        accumulators[1].setBennetParam(refPrefCenter,span);
    }
    
    public void setMeter(MeterVirial newMeter, int iBox) {
        if (accumulators[iBox] != null) {
            // we need a new accumulator so nuke the old one now.
            if (accumulatorPumps[iBox] != null) {
                integrators[iBox].getEventManager().removeListener(accumulatorPumps[iBox]);
                accumulatorPumps[iBox] = null;
            }
            accumulators[iBox] = null;
        }
        meters[iBox] = newMeter;
        newMeter.setBox(box[iBox]);
    }

    public void setAccumulator(AccumulatorVirialOverlapSingleAverage newAccumulator, int iBox) {
        accumulators[iBox] = newAccumulator;
        accumulators[iBox].setBlockSize(blockSize);
        if (accumulatorPumps[iBox] == null) {
            accumulatorPumps[iBox] = new DataPumpListener(meters[iBox],newAccumulator);
            integrators[iBox].getEventManager().addListener(accumulatorPumps[iBox]);
        }
        else {
            accumulatorPumps[iBox].setDataSink(newAccumulator);
        }
        if (integratorOS != null) {
            dsvo = new DataSourceVirialOverlap(accumulators[0],accumulators[1]);
            integratorOS.setReferenceFracSource(dsvo);
        }
    }
    
    public void setAccumulatorBlockSize(long newBlockSize) {
        blockSize = newBlockSize;
        for (int i=0; i<2; i++) {
            accumulators[i].setBlockSize(newBlockSize);
        }
        // reset the integrator so that it will re-adjust step frequency
        // and ensure it will take enough data for both ref and target
        integratorOS.reset();
    }

    /**
     * Sets the number of alpha values used for the production stage of the
     * simulation
     */
    public void setNumAlpha(int newNumAlpha) {
        numAlpha = newNumAlpha;
    }

    /**
     * Returns the number of alpha values used for the production stage of the
     * simulation
     */
    public int getNumAlpha() {
        return numAlpha;
    }

    public void setRefPref(double newRefPref) {
        System.out.println("setting ref pref (explicitly) to "+newRefPref);
        setAccumulator(new AccumulatorVirialOverlapSingleAverage(1,true),0);
        setAccumulator(new AccumulatorVirialOverlapSingleAverage(1,false),1);
        setRefPref(newRefPref,1);
    }
    
    public void initRefPref(String fileName, long initSteps) {
        // use the old refpref value as a starting point so that an initial
        // guess can be provided
        double oldRefPref = refPref;
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
                setAccumulator(new AccumulatorVirialOverlapSingleAverage(numAlpha,true),0);
                setAccumulator(new AccumulatorVirialOverlapSingleAverage(numAlpha,false),1);
                setRefPref(refPref,1);
            }
            catch (IOException e) {
                // file not there, which is ok.
            }
        }
        
        if (refPref == -1) {
            for (int i=0; i<2; i++) {
                integrators[i].getMoveManager().setEquilibrating(true);
            }

            long oldBlockSize = blockSize;
            // 1000 blocks
            long newBlockSize = initSteps*integratorOS.getNumSubSteps()/1000;
            if (newBlockSize < 1000) {
                // make block size at least 1000, even if it means fewer blocks
                newBlockSize = 1000;
            }
            if (newBlockSize > 1000000) {
                // needs to be an int.  1e6 steps/block is a bit crazy.
                newBlockSize = 1000000;
            }
            setAccumulatorBlockSize(newBlockSize);
            setAccumulator(new AccumulatorVirialOverlapSingleAverage(21,true),0);
            setAccumulator(new AccumulatorVirialOverlapSingleAverage(21,false),1);
            setRefPref(oldRefPref,30);
            boolean adjustable = integratorOS.isAdjustStepFraction();
            if (adjustable) {
                // we do this initialization to
                // 1. find alpha
                // 2. get molecules out of their starting configuration
                // 3. find optimal mc move step sizes
                // all of these are about as hard in the reference as in the target system
                // so force integratorOS to run both systems equally.
                integratorOS.setRefStepFraction(0.5);
                integratorOS.setAdjustStepFraction(false);
            }
            ai.setMaxSteps(initSteps);
            ai.actionPerformed();
            if (adjustable) {
                integratorOS.setAdjustStepFraction(true);
            }

            int newMinDiffLoc = dsvo.minDiffLocation();
            refPref = accumulators[0].getBennetAverage(newMinDiffLoc)
                /accumulators[1].getBennetAverage(newMinDiffLoc);
            System.out.println("setting ref pref to "+refPref);
            if (Double.isInfinite(refPref) || Double.isNaN(refPref)) {
                throw new RuntimeException("oops "+accumulators[0].getBennetAverage(newMinDiffLoc)+" "+accumulators[1].getBennetAverage(newMinDiffLoc));
            }
            setAccumulatorBlockSize(oldBlockSize);
            setAccumulator(new AccumulatorVirialOverlapSingleAverage(15,true),0);
            setAccumulator(new AccumulatorVirialOverlapSingleAverage(15,false),1);
            setRefPref(refPref,4);
            for (int i=0; i<2; i++) {
                integrators[i].reset();
            }
            // set refPref back to -1 so that later on we know that we've been looking for
            // the appropriate value
            refPref = -1;
        }

    }
    
    public void equilibrate(String fileName, long initSteps) {
        // run a short simulation to get reasonable MC Move step sizes and
        // (if needed) narrow in on a reference preference
        ai.setMaxSteps(initSteps);
        long oldBlockSize = blockSize;
        // 1000 blocks
        long newBlockSize = initSteps*integratorOS.getNumSubSteps()/1000;
        if (newBlockSize < 1000) {
            // make block size at least 1000, even if it means fewer blocks
            newBlockSize = 1000;
        }
        if (newBlockSize > 1000000) {
            // needs to be an int.  1e6 steps/block is a bit crazy.
            newBlockSize = 1000000;
        }
        setAccumulatorBlockSize((int)newBlockSize);
        for (int i=0; i<2; i++) {
            integrators[i].getMoveManager().setEquilibrating(true);
        }
        boolean adjustable = integratorOS.isAdjustStepFraction();
        if (adjustable) {
            // we do this initialization to
            // 1. find alpha
            // 2. get molecules out of their starting configuration
            // 3. find optimal mc move step sizes
            // all of these are about as hard in the reference as in the target system
            // so force integratorOS to run both systems equally.
            integratorOS.setRefStepFraction(0.5);
            integratorOS.setAdjustStepFraction(false);
        }
        ai.actionPerformed();
        if (adjustable) {
            integratorOS.setAdjustStepFraction(true);
        }

        if (refPref == -1) {
            int newMinDiffLoc = dsvo.minDiffLocation();
            refPref = accumulators[0].getBennetAverage(newMinDiffLoc)
                /accumulators[1].getBennetAverage(newMinDiffLoc);
            System.out.println("setting ref pref to "+refPref+" ("+newMinDiffLoc+")");
            if (Double.isInfinite(refPref) || Double.isNaN(refPref)) {
                throw new RuntimeException("oops "+accumulators[0].getBennetAverage(newMinDiffLoc)+" "+accumulators[1].getBennetAverage(newMinDiffLoc));
            }
            setAccumulator(new AccumulatorVirialOverlapSingleAverage(numAlpha,true),0);
            setAccumulator(new AccumulatorVirialOverlapSingleAverage(numAlpha,false),1);
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
        setAccumulatorBlockSize(oldBlockSize);
        for (int i=0; i<2; i++) {
            integrators[i].getMoveManager().setEquilibrating(false);
        }
    }
    
    private static final long serialVersionUID = 1L;
	public DataSourceVirialOverlap dsvo;
    public AccumulatorVirialOverlapSingleAverage[] accumulators;
    protected DataPumpListener[] accumulatorPumps;
	protected final ClusterWeight[] sampleClusters;
    public BoxCluster[] box;
    public IntegratorMC[] integrators;
    public MCMoveBoxStep[] mcMoveRotate;
    public MCMoveBoxStep[] mcMoveTranslate;
    public MCMoveBoxStep[] mcMoveWiggle;
    public MeterVirial[] meters;
    public ActivityIntegrate ai;
    public IntegratorOverlap integratorOS;
    public double refPref;
    protected long blockSize;
    protected int numAlpha = 1;
}
