/*
 * History
 * Created on Aug 11, 2004 by kofke
 */
package etomica;

/**
 * @author kofke
 *
 * Interface for an AtomIterator that can be conditioned with
 * target and basis atoms.
 */
public interface AtomsetIteratorBasisDependent extends AtomsetIteratorTargetable {

	/**
	 * Identifies the atoms that form the basis for iteration, such that
	 * the childList atoms of those given will form the iterates.
	 * @param atoms The basis atoms; a null or zero-length array will
	 * condition the iterator to give no iterates until a valid basis
	 * is specified via another call to this method.
	 */
    public void setBasis(Atom[] atoms);
    
    /**
     * Indicates the size of the basis needed to set the iterator.
     * The length of the array given to setBasis should be this value.
     * @return the size of the basis for this iterator.
     */
    public int basisSize();

}
