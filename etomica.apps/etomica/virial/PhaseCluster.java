package etomica.virial;

import etomica.atom.AtomArrayList;
import etomica.atom.AtomLeaf;
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
    public void trialNotify() {
        // atom(s) have been moved.  leave cPairSet as is and update
        // cPairTrialSet and set a flag to use it.
        isTrial = true;
        // increase ID to notify clusters to recalculate value
        if(cPairSet == null) {
            // assume 1 species
            AtomArrayList molecules = ((AtomTreeNodeGroup)((AtomTreeNodeGroup)getSpeciesMaster().getNode()).childList.get(0).getNode()).childList;
            if (molecules.get(0) instanceof AtomLeaf) {
                cPairSet = new CoordinatePairLeafSet(molecules,space);
                cPairTrialSet = new CoordinatePairLeafSet(molecules,space);
            }
            else {
                cPairSet = new CoordinatePairMoleculeSet(molecules,space);
                cPairTrialSet = new CoordinatePairMoleculeSet(molecules,space);
            }
            aPairSet = new AtomPairSet(molecules);
        }

        cPairTrialSet.reset();
    }
	
    /**
     * Informs the phase that the trial was accepted so it will keep the new 
     * coordinate pairs.
     */
	public void acceptNotify() {
        if (!isTrial) {
            throw new IllegalStateException("you weren't in a trial!");
        }
	    // move was accepted.  swap out trial cPairSet and cPairTrialSet since
		// cPairTrialSet is already up-to-date
		isTrial = false;
		cPairSetTmp = cPairSet;
		cPairSet = cPairTrialSet;
		cPairTrialSet = cPairSetTmp;
	}
	
    /**
     * Informs the phase that the trial was rejected so it will go back to 
     * the old coordinate pairs.
     */
	public void rejectNotify() {
        if (!isTrial) {
            throw new IllegalStateException("you weren't in a trial!");
        }
		// move was rejected.  stop using cPairTrialSet.
		isTrial = false;
	}
	
    private static final long serialVersionUID = 1L;
	private boolean isTrial;
	private CoordinatePairSet cPairSet, cPairTrialSet, cPairSetTmp;
    private AtomPairSet aPairSet;
	private final ClusterWeight sampleCluster;
}