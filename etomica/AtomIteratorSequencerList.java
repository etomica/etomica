/*
 * History
 * Created on Aug 26, 2004 by kofke
 */
package etomica;

/**
 * Iterator for looping through the sequence list relative to a
 * specified atom.
 */
public class AtomIteratorSequencerList extends AtomIteratorAdapter implements AtomIteratorAtomDependent {

	public AtomIteratorSequencerList() {
		super(new AtomIteratorList());
		listIterator = (AtomIteratorList)iterator;
		listIterator.unset();
	}
	
	public void reset() {
		if(firstAtom == null) return;
		listIterator.reset();
		if(skippingFirst) listIterator.next();
	}

	/* (non-Javadoc)
	 * @see etomica.AtomIterator#allAtoms(etomica.AtomActive)
	 */
	public void allAtoms(AtomActive action) {
		//TODO this method
		super.allAtoms(action);
	}
	/* (non-Javadoc)
	 * @see etomica.AtomIterator#contains(etomica.Atom)
	 */
	public boolean contains(Atom atom) {
		return skippingFirst ? (iterator.contains(atom) && (atom == firstAtom)) 
							 : iterator.contains(atom);
	}
	/* (non-Javadoc)
	 * @see etomica.AtomIterator#size()
	 */
	public int size() {
		return skippingFirst ? listIterator.size() : listIterator.size() - 1;
	}
	/* (non-Javadoc)
	 * @see etomica.AtomIteratorPhaseDependent#setPhase(etomica.Phase)
	 */
	public void setAtom(Atom atom) {
		firstAtom = atom;
		if(atom != null) listIterator.setFirst(atom.seq);
		else listIterator.unset();
	}

	/**
	 * @return Returns the skippingFirst.
	 */
	public boolean isSkippingFirst() {
		return skippingFirst;
	}
	/**
	 * @param skippingFirst The skippingFirst to set.
	 */
	public void setSkippingFirst(boolean skippingFirst) {
		this.skippingFirst = skippingFirst;
	}
	
	public void setIterationDirection(IteratorDirective.Direction direction) {
		listIterator.setIterationDirection(direction);
	}
	public IteratorDirective.Direction getIterationDirection() {
		return listIterator.getIterationDirection();
	}
	
	private final AtomIteratorList listIterator;
	private boolean skippingFirst = false;
	private Atom firstAtom = null;

}
