/*
 * History
 * Created on Aug 11, 2004 by kofke
 */
package etomica;

/**
 * @author kofke
 *
 * Interface for an AtomIterator that can be reset with an IteratorDirective.
 */
public interface AtomIteratorDirectable extends AtomIterator {

    public void reset(IteratorDirective id);

}
