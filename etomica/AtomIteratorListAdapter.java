/*
 * History
 * Created on Aug 13, 2004 by kofke
 */
package etomica;

/**
 * @author kofke
 * Adapater class that wraps an AtomIteratorList to implement the
 * methods of an iterator.  Subclasses must define the setPhase method
 * so that it configures in an appropriate way for a given phase.
 */
public abstract class AtomIteratorListAdapter implements AIPhaseDependent {

	protected final AIAtomListDependent listIterator;
	
	/**
	 * 
	 */
	public AtomIteratorListAdapter(AIAtomListDependent aiali) {
		this.listIterator = aiali;
	}

	/* (non-Javadoc)
	 * @see etomica.AtomIterator#contains(etomica.Atom)
	 */
	public final boolean contains(Atom atom) {
		return listIterator.contains(atom);
	}

	/* (non-Javadoc)
	 * @see etomica.AtomIterator#hasNext()
	 */
	public final boolean hasNext() {
		return listIterator.hasNext();
	}

	/* (non-Javadoc)
	 * @see etomica.AtomIterator#reset()
	 */
	public final Atom reset() {
		return listIterator.reset();
	}

	/* (non-Javadoc)
	 * @see etomica.AtomIterator#unset()
	 */
	public final void unset() {
		listIterator.unset();
	}

	/* (non-Javadoc)
	 * @see etomica.AtomIterator#next()
	 */
	public final Atom next() {
		return listIterator.next();
	}

	/* (non-Javadoc)
	 * @see etomica.AtomIterator#all(etomica.Atom, etomica.IteratorDirective, etomica.AtomActive)
	 */
	public final void all(Atom basis, IteratorDirective id, AtomActive action) {
		listIterator.all(basis, id, action);
	}

	/* (non-Javadoc)
	 * @see etomica.AtomIterator#size()
	 */
	public final int size() {
		return listIterator.size();
	}
	
	public abstract void setPhase(Phase phase);
	
	
}
