/*
 * History
 * Created on Oct 19, 2004 by kofke
 */
package etomica;

/**
 * Interface for an iterator that requires specification of an 
 * AtomList for its operation.
 */
public interface AtomIteratorListDependent extends AtomIterator {

	/**
	 * Sets the list for iteration.
	 */
	public abstract void setList(AtomList list);
	
	/**
	 * Returns the list currently used for iteration.
	 */
	public AtomList getList();
}
