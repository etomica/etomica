package etomica;

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
	public PotentialN(SimulationElement parent) {
		super(Integer.MAX_VALUE, parent);
	}

	/**
	 * @see etomica.Potential#calculate(etomica.AtomSet, etomica.IteratorDirective, etomica.PotentialCalculation)
	 */
//	public void calculate(AtomSet basis, IteratorDirective id, PotentialCalculation pc) {
//		//still abstract
//			if(!enabled) return;
//			pc.set(this).actionPerformed(phase);
//	}

	public abstract double energy(AtomSet atomSet);

	/**
	 * Sets the iterator for this potential, which in most cases does not
	 * require any iterator.  Default is AtomSetIterator.NULL.
	 * @see etomica.Potential#setIterator(AtomSetIterator)
	 */
	public void setIterator(AtomSetIterator iterator) {
		this.iterator = iterator;
	}
    
	/**
	 * Returns the iterator last defined via the setIterator method.  Default is
	 * AtomSetIterator.NULL if none was previously set.
	 * @see etomica.Potential#getIterator()
	 */
	public AtomSetIterator getIterator() {
		return iterator;
	}

	private AtomSetIterator iterator = AtomSetIterator.NULL;

}
