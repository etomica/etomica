package etomica.atom.iterator;

import etomica.atom.IAtom;

/**
 * Interface for an iterator that can be targeted at one or more atoms.
 * Setting a target causes the iterator to produce only those of its iterates 
 * that contain the target atom(s).  A typical use would be to set an atom
 * as a target for a pair iterator, so the iterator would form only those pairs
 * that contain the targeted atom.  The argument can be null, indicating no 
 * target.
 */
public interface AtomsetIteratorTargetable extends AtomsetIterator {
	public void setTarget(IAtom targetAtom);
}