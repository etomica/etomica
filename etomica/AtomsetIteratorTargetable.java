/*
 * History
 * Created on Aug 31, 2004 by kofke
 */
package etomica;

/**
 * Interface for an iterator that can be targeted at one or more atoms.
 * Setting a target causes the iterator to produce only those of its iterates 
 * that contain the target atom(s).  A typical use would be to set an atom
 * as a target for a pair iterator, so the iterator would form only those pairs
 * that contain the targeted atom.  A null or zero-length array argument can
 * be given to indicate that that the iterator should gives its full complement
 * of iterates.
 */
public interface AtomsetIteratorTargetable extends AtomsetIterator {
	public void setTarget(Atom[] targetAtoms);
}