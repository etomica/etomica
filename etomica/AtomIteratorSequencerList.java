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

	/**
	 * Constructs new class with hasNext as false.  Must
	 * invoke setAtom and reset before beginning iteration.
	 */
	public AtomIteratorSequencerList() {
		super(new AtomIteratorList());
		listIterator = (AtomIteratorList)iterator;
		listIterator.unset();
	}
	
	/**
	 * Overrides superclass reset to ensure no reset is performed
	 * if a firstAtom has not been identified.  Otherwise readies
	 * for iteration beginning with firstAtom, or the one following
	 * it if skippingFirst is true.
	 */
	public void reset() {
		if(firstAtom == null) return;
		listIterator.reset();
		if(skippingFirst) listIterator.next();
	}

	/* (non-Javadoc)
	 * @see etomica.AtomIterator#allAtoms(etomica.AtomActive)
	 */
	public void allAtoms(AtomsetActive action) {
		if(skippingFirst) listIterator.allAtoms(skipFirstWrapper(action));
		else listIterator.allAtoms(action);
	}
	/* (non-Javadoc)
	 * @see etomica.AtomIterator#contains(etomica.Atom)
	 */
	public boolean contains(Atom[] atom) {
		return skippingFirst ? (atom[0] != firstAtom && listIterator.contains(atom)) 
							 : listIterator.contains(atom);
	}
	/* (non-Javadoc)
	 * @see etomica.AtomIterator#size()
	 */
	public int size() {
		return skippingFirst ? listIterator.size() : listIterator.size() - 1;
	}
	/* (non-Javadoc)
	 * @see etomica.AtomsetIteratorPhaseDependent#setPhase(etomica.Phase)
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

	/**
	 * Constructs an AtomActive class that wraps another, with the
	 * effect of preventing it from performing its action on firstAtom.
	 * Used by the allAtoms method if skippingFirst.
	 */
	private AtomsetActive skipFirstWrapper(final AtomsetActive realAction) {
        return new AtomsetActive() {
                //fake action for first iterate
                private AtomsetActive action = new AtomsetActive() {
                        public void actionPerformed(Atom[] atom) {action = realAction;}
                };
                public void actionPerformed(Atom[] atom) {
                        action.actionPerformed(atom);
                }
        };
}

}
