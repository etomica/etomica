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

}