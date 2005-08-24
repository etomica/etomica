package etomica.atom;


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
     * Returns true for null atom also.
     */
    public static final AtomFilter ACCEPT_ALL = new AtomFilter() {
        public boolean accept(Atom a) {return true;}
    };
    
    /**
     * Static instance of a filter that rejects all atoms.
     * Returns false for null atom also.
     */
    public static final AtomFilter ACCEPT_NONE = new AtomFilter() {
        public boolean accept(Atom a) {return false;}
    };
    
}