/*
 * History
 * Created on Aug 13, 2004 by kofke
 */
package etomica;

/**
 * Adapater class that wraps another AtomIterator to implement the
 * methods of an iterator, while adding a PhaseDependent interface.  
 * Subclasses must define the setPhase method so that it configures 
 * the wrapped iterator to perform the desired iteration within
 * the given phase.
 */
public abstract class AtomIteratorAdapter implements AIPhaseDependent {

	protected final AtomIterator iterator;
	
	/**
	 * 
	 */
	public AtomIteratorAdapter(AtomIterator iterator) {
		this.iterator = iterator;
	}

	/* (non-Javadoc)
	 * @see etomica.AtomIterator#contains(etomica.Atom)
	 */
	public final boolean contains(Atom atom) {
		return iterator.contains(atom);
	}

	/* (non-Javadoc)
	 * @see etomica.AtomIterator#hasNext()
	 */
	public final boolean hasNext() {
		return iterator.hasNext();
	}

	/* (non-Javadoc)
	 * @see etomica.AtomIterator#reset()
	 */
	public final void reset() {
		iterator.reset();
	}

	/* (non-Javadoc)
	 * @see etomica.AtomIterator#unset()
	 */
	public final void unset() {
		iterator.unset();
	}

	/* (non-Javadoc)
	 * @see etomica.AtomIterator#next()
	 */
	public final Atom next() {
		return iterator.next();
	}
	
	public final Atom peek() {
		return iterator.peek();
	}

	/* (non-Javadoc)
	 * @see etomica.AtomIterator#all(etomica.Atom, etomica.IteratorDirective, etomica.AtomActive)
	 */
	public final void allAtoms(AtomActive action) {
		iterator.allAtoms(action);
	}

	/* (non-Javadoc)
	 * @see etomica.AtomIterator#size()
	 */
	public final int size() {
		return iterator.size();
	}
	
	public abstract void setPhase(Phase phase);
	
}
