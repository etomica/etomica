/*
 * History
 * Created on Aug 11, 2004 by kofke
 */
package etomica;

/**
 * @author kofke
 *
 * Interface for an AtomIterator that can be conditioned with an
 * IteratorDirective and a set of basis atoms.
 */
public interface AtomsetIteratorDirectable extends AtomsetIterator {

    public void setDirective(IteratorDirective id);
    
    public void setBasis(Atom[] atoms);
    
    /**
     * Indicates the size of the basis needed to set the iterator.
     * The length of the array given to setBasis should be this value.
     * @return the size of the basis for this iterator.
     */
    public int basisSize();

}
