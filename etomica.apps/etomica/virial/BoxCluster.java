package etomica.virial;

import etomica.atom.AtomSet;
import etomica.box.Box;
import etomica.simulation.ISimulation;
import etomica.space.BoundaryRectangularNonperiodic;

/**
 * @author kofke
 *
 * Extension of Box that forms and holds a PairSet instance for all of the
 * atoms in the box.  Also instantiates box with a NONE boundary type.
 */
public class BoxCluster extends Box {

	/**
	 * Constructor for BoxCluster.
	 */
	public BoxCluster(ISimulation sim, ClusterWeight cluster) {
		super(new BoundaryRectangularNonperiodic(sim.getSpace(), sim.getRandom()));
        sampleCluster = cluster;
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
     * returns the cluster used for sampling in this box
	 */
    public ClusterWeight getSampleCluster() {
        return sampleCluster;
    }

    /**
     * Inform the box that a trial move has been made so it can update
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
            AtomSet molecules = getMoleculeList();
            AtomSet leafAtoms = getLeafList();
            if (molecules.getAtomCount() == leafAtoms.getAtomCount()) {
                cPairSet = new CoordinatePairLeafSet(leafAtoms,space);
                cPairTrialSet = new CoordinatePairLeafSet(leafAtoms,space);
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
     * Informs the box that the trial was accepted so it will keep the new 
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
     * Informs the box that the trial was rejected so it will go back to 
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