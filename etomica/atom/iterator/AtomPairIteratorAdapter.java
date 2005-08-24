/*
 * History
 * Created on Aug 13, 2004 by kofke
 */
package etomica.atom.iterator;

import etomica.atom.AtomPair;
import etomica.atom.AtomSet;

/**
 * Adapater class that wraps another AtomPairIterator to implement the
 * methods of an iterator.  In the default case, all methods of the
 * adapter merely front the methods of the wrapped iterator.
 * Typically, subclasses will introduce additional
 * methods that configure the wrapped iterator in a particular way to
 * prepare it for iteration.  Subclasses may also override individual
 * methods as appropriate to their implementation.
 */
public abstract class AtomPairIteratorAdapter extends AtomsetIteratorAdapter implements AtomPairIterator {

	protected final AtomPairIterator atomPairIterator;
	
	/**
	 * Constructor takes the wrapped iterator as an argument. The
	 * iterator is final and cannot be subsequently replaced by
	 * another iterator.
	 */
	public AtomPairIteratorAdapter(AtomPairIterator iterator) {
		super(iterator);
		atomPairIterator = iterator;
	}
	
    public final AtomSet next() {
        return nextPair();
    }
    
	public AtomPair nextPair() {
		return atomPairIterator.nextPair();
	}
}
