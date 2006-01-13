package etomica.virial;

import etomica.atom.Atom;
import etomica.atom.AtomArrayList;
import etomica.atom.AtomTreeNodeGroup;
import etomica.phase.Phase;
import etomica.simulation.Simulation;
import etomica.space.BoundaryRectangularNonperiodic;

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
        setBoundary(new BoundaryRectangularNonperiodic(sim.space));
	}
	
    /**
     * returns the current coordinate pair set
     */
	public CoordinatePairSet getCPairSet() {
		return isTrial ? cPairTrialSet : cPairSet;
	}
    
    public AtomPairSet getAPairSet() {
        return aPairSet;
    }
    
	/**
     * returns the cluster used for sampling in this phase
	 */
    public ClusterWeight getSampleCluster() {
        return sampleCluster;
    }

    /**
     * Inform the phase that a trial move has been made so it can update
     * the coordinate pairs.  If molecule is not null, only coordinate pairs 
     * containing that atom are updated.
     */
    public void trialNotify(Atom molecule) {
        // atom(s) have been moved.  leave cPairSet as is and update
        // cPairTrialSet and set a flag to use it.
        isTrial = true;
        // increase ID to notify clusters to recalculate value
        if(cPairSet == null) {
            // assume 1 species
            AtomArrayList molecules = ((AtomTreeNodeGroup)((AtomTreeNodeGroup)getSpeciesMaster().node).childList.get(0).node).childList;
            cPairSet = new CoordinatePairSet(molecules,space);
            cPairTrialSet = new CoordinatePairSet(molecules,space);
            aPairSet = new AtomPairSet(molecules);
        }

        if (molecule != null) {
            int dirtyAtom = molecule.node.getIndex();
            // copy old cPairSet to trial cPairSet
            cPairTrialSet.E(cPairSet);
            cPairTrialSet.dirtyAtom = dirtyAtom;
        }
        else {
            cPairTrialSet.dirtyAtom = -1;
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
    private AtomPairSet aPairSet;
	private final ClusterWeight sampleCluster;
}