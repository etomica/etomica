/*
 * History
 * Created on Aug 13, 2004 by kofke
 */
package etomica;

/**
 * Adapater class that wraps another AtomsetIterator to implement the
 * methods of an iterator.  In the default case, all methods of the
 * adapter merely front the methods of the wrapped iterator.
 * Typically, subclasses will introduce additional
 * methods that configure the wrapped iterator in a particular way to
 * prepare it for iteration.  Subclasses may also override individual
 * methods as appropriate to their implementation.
 */
public abstract class AtomsetIteratorAdapter implements AtomsetIterator {

	protected final AtomsetIterator iterator;
	
	/**
	 * Constructor takes the wrapped iterator as an argument. The
	 * iterator is final and cannot be subsequently replaced by
	 * another iterator.
	 */
	public AtomsetIteratorAdapter(AtomsetIterator iterator) {
		this.iterator = iterator;
	}

	/* (non-Javadoc)
	 * @see etomica.AtomIterator#contains(etomica.Atom)
	 */
	public boolean contains(Atom[] atom) {
		return iterator.contains(atom);
	}

	/* (non-Javadoc)
	 * @see etomica.AtomIterator#hasNext()
	 */
	public boolean hasNext() {
		return iterator.hasNext();
	}

	/* (non-Javadoc)
	 * @see etomica.AtomIterator#reset()
	 */
	public void reset() {
		iterator.reset();
	}

	/* (non-Javadoc)
	 * @see etomica.AtomIterator#unset()
	 */
	public void unset() {
		iterator.unset();
	}

	/* (non-Javadoc)
	 * @see etomica.AtomIterator#next()
	 */
	public Atom[] next() {
		return iterator.next();
	}
	
	public Atom[] peek() {
		return iterator.peek();
	}
	
	/* (non-Javadoc)
	 * @see etomica.AtomIterator#all(etomica.Atom, etomica.IteratorDirective, etomica.AtomActive)
	 */
	public void allAtoms(AtomsetActive action) {
		iterator.allAtoms(action);
	}

	/* (non-Javadoc)
	 * @see etomica.AtomIterator#size()
	 */
	public int size() {
		return iterator.size();
	}
	
	public final int nBody() {
		return iterator.nBody();
	}
}
