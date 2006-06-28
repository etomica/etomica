package etomica.potential;

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
	public PotentialN(Space space){
		this(space, Integer.MAX_VALUE);
	}
	
	public PotentialN(Space space, int nBody) {
		super(nBody, space);
	}
}
