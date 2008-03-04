/*
 * History
 * Created on Aug 13, 2004 by kofke
 */
package etomica.atom.iterator;

import etomica.action.AtomsetAction;
import etomica.api.IAtomSet;

/**
 * Adapater class that wraps another AtomsetIterator to implement the
 * methods of an iterator.  In the default case, all methods of the
 * adapter merely front the methods of the wrapped iterator.
 * Typically, subclasses will introduce additional
 * methods that configure the wrapped iterator in a particular way to
 * prepare it for iteration.  Subclasses may also override individual
 * methods as appropriate to their implementation.
 */
public abstract class AtomsetIteratorAdapter implements AtomsetIterator, java.io.Serializable {

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
	public IAtomSet next() {
		return iterator.next();
	}
	
	/* (non-Javadoc)
	 * @see etomica.AtomIterator#all(etomica.Atom, etomica.IteratorDirective, etomica.AtomActive)
	 */
	public void allAtoms(AtomsetAction action) {
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
