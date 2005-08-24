/*
 * History
 * Created on Dec 3, 2004 by kofke
 */
package etomica.atom.iterator;

import etomica.atom.AtomIterator;
import etomica.atom.AtomPairIterator;

/**
 * Interface for an atom pair iterator that is formed from two atom iterators.
 */
public interface ApiComposite extends AtomPairIterator {

    /**
     * @return inner-loop iterator
     */
    public AtomIterator getInnerIterator();
    
    /**
     * @return outer-loop iterator
     */
    public AtomIterator getOuterIterator();
}
