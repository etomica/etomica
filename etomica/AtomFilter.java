package etomica;

/**
 * Interface for a class that screens atoms according
 * to some criterion.
 */
 
 /* History of changes
  * 09/07/02 (DAK) new
  */
  
public interface AtomFilter {
    
    /**
     * Returns true if atom is passes test of filter.
     */
    public boolean accept(Atom a);

//------------ end of interface ---------------//


    /**
     * Static instance of a filter that accepts all atoms.
     */
    public static final AtomFilter ALL = new All();
    
    /**
     * Static instance of a filter that rejects all atoms.
     */
    public static final AtomFilter NONE = new None();
    
    /**
     * Filter that returns true for all atoms.
     */
    static class All implements AtomFilter {
        public boolean accept(Atom a) {return true;}
    };
    
    /**
     * Filter that returns false for all atoms.
     */
    static class None implements AtomFilter {
        public boolean accept(Atom a) {return false;}
    }
    
}