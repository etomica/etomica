package etomica.atom;

import etomica.AtomPair;

/**
 * Interface for a class that screens atom pairs according
 * to some criterion.
 */
 
public interface AtomPairFilter {
    
    /**
     * Returns true if atom pair passes test of filter.
     */
    public boolean accept(AtomPair pair);

//------------ end of interface ---------------//


    /**
     * Static instance of a filter that accepts all pairs.
     * Returns true for null pair also.
     */
    public static final AtomPairFilter ACCEPT_ALL = new AtomPairFilter() {
        public boolean accept(AtomPair pair) {return true;}
    };
    
    /**
     * Static instance of a filter that rejects all pairs.
     * Returns false for null pair also.
     */
    public static final AtomPairFilter ACCEPT_NONE = new AtomPairFilter() {
        public boolean accept(AtomPair pair) {return false;}
    };
    
}