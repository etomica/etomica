package etomica.potential;

import etomica.atom.AtomSet;
import etomica.nbr.NeighborCriterion;
import etomica.phase.Phase;
import etomica.space.Space;

/**
 * @author kofke
 *
 * General potential that depends on positions of all N molecules, or is
 * otherwise not naturally expressed as a single-, pair-, etc-body potential.
 */

/* History
 * 08/29/03 (DAK) new; introduced for etomica.research.nonequilwork.PotentialOSInsert
 */
public abstract class PotentialN extends Potential {

	/**
	 * Constructor for PotentialN.
	 * @param sim
	 */
	public PotentialN(int nBody, Space space) {
		super(nBody, space);
	}

    public void setCriterion(NeighborCriterion criterion) {
        this.criterion = criterion;
    }

    public NeighborCriterion getCriterion() {
        return criterion;
    }

    protected NeighborCriterion criterion = NeighborCriterion.ALL;
}
