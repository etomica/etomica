package etomica.atom;


/**
 * Interface for a class that screens atoms according
 * to some criterion.
 */
public interface AtomFilter {
    
    /**
     * Returns true if atom is passes test of filter.
     */
    public boolean accept(IAtom a);

}