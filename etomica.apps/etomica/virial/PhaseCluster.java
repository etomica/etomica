package etomica.virial;

import etomica.Phase;
import etomica.Simulation;
import etomica.space.BoundaryNone;

/**
 * @author kofke
 *
 * Extension of Phase that forms and holds a PairSet instance for all of the
 * atoms in the phase.  Also instantiates phase with a NONE boundary type.
 */
public class PhaseCluster extends Phase {

	/**
	 * Constructor for PhaseCluster.
	 */
	public PhaseCluster(Simulation sim, ClusterWeight cluster) {
		super(sim);
        sampleCluster = cluster;
        setBoundary(new BoundaryNone(sim.space));
        // this is a bit contorted, but config needs to know the phase before 
        // the coordinates are initialized
        ConfigurationCluster config = new ConfigurationCluster(sim.space);
        config.setPhase(this);
        setConfiguration(config);
	}
	
    /**
     * returns the current coordinate pair set
     */
	public CoordinatePairSet getCPairSet() {
		return isTrial ? cPairTrialSet : cPairSet;
	}
    
	/**
     * returns the cluster used for sampling in this phase
	 */
    public ClusterWeight getSampleCluster() {
        return sampleCluster;
    }

    /**
     * Inform the phase that a trial move has been made so it can update
     * the coordinate pairs.
     */
	public void trialNotify() {
		// atom(s) have been moved.  leave cPairSet as is and update
		// cPairTrialSet and set a flag to use it.
		isTrial = true;
		// increase ID to notify clusters to recalculate value
        if(cPairSet == null) {
            cPairSet = new CoordinatePairSet(speciesMaster.atomList,space);
            cPairTrialSet = new CoordinatePairSet(speciesMaster.atomList,space);
        }
		cPairTrialSet.reset();
	}
	
    /**
     * Informs the phase that the trial was accepted so it will keep the new 
     * coordinate pairs.
     */
	public void acceptNotify() {
	    // move was accepted.  swap out trial cPairSet and cPairTrialSet since
		// cPairTrialSet is already up-to-date
		isTrial = false;
		cPairSetTmp = cPairSet;
		cPairSet = cPairTrialSet;
		cPairTrialSet = cPairSetTmp;
	}
	
    /**
     * Informs the phase that the trial was accepted so it will go back to 
     * the old coordinate pairs.
     */
	public void rejectNotify() {
		// move was rejected.  stop using cPairTrialSet.
		isTrial = false;
	}
	
	private boolean isTrial;
	private CoordinatePairSet cPairSet, cPairTrialSet, cPairSetTmp;
	private final ClusterWeight sampleCluster;
}