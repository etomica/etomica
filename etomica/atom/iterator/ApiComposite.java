package etomica.atom.iterator;


/**
 * Interface for an atom pair iterator that is formed from two atom iterators.
 */
public interface ApiComposite extends AtomsetIterator {

    /**
     * @return inner-loop iterator
     */
    public AtomIterator getInnerIterator();
    
    /**
     * @return outer-loop iterator
     */
    public AtomIterator getOuterIterator();
}
