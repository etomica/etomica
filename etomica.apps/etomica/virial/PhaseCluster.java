package etomica.virial;

import etomica.Phase;
import etomica.Space;
import etomica.space.BoundaryNone;
import etomica.space3d.Space3D;

/**
 * @author kofke
 *
 * Extension of Phase that forms and holds a PairSet instance for all of the
 * atoms in the phase.  Also instantiates phase with a NONE boundary type.
 */
public class PhaseCluster extends Phase {

	/**
	 * Constructor for PhaseCluster.
	 * @param parent
	 */
	public PhaseCluster(Space space, ClusterWeight cluster) {
		super(space);
        sampleCluster = cluster;
        setBoundary(new BoundaryNone(space));
	}
	
	public CoordinatePairSet getCPairSet() {
		if(cPairSet == null && speciesMaster.atomList.size() > 0) {
			cPairSet = new CoordinatePairSet(speciesMaster.atomList,new Space3D());
			cPairTrialSet = new CoordinatePairSet(speciesMaster.atomList,new Space3D());
			// coordinate pairs are "reset" in the constructor, so we don't need to do it here
		}
		return isTrial ? cPairTrialSet : cPairSet;
	}
    public ClusterWeight getSampleCluster() {
        return sampleCluster;
    }
	
	public void trialNotify() {
		// atom(s) have been moved.  leave cPairSet as is and update
		// cPairTrialSet and set a flag to use it.
		isTrial = true;
		// increase ID to notify clusters to recalculate value
		cPairTrialSet.reset();
	}
	
	public void acceptNotify() {
	    // move was accepted.  swap out trial cPairSet and cPairTrialSet since
		// cPairTrialSet is already up-to-date
		isTrial = false;
		cPairSetTmp = cPairSet;
		cPairSet = cPairTrialSet;
		cPairTrialSet = cPairSetTmp;
	}
	
	public void rejectNotify() {
		// move was rejected.  stop using cPairTrialSet.
		isTrial = false;
	}
	
	private boolean isTrial;
	private CoordinatePairSet cPairSet, cPairTrialSet, cPairSetTmp;
	private final ClusterWeight sampleCluster;
}
