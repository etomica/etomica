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

}
