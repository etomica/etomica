package etomica;

/**
 * Methods for properties obtained for a soft, differentiable pair potential.
 *
 * @author David Kofke
 */

public abstract class Potential2Soft extends Potential2 {
    
    
    public Potential2Soft(PotentialGroup parent) {
        super(parent);
    }
    public Potential2Soft(PotentialGroup parent, PotentialTruncation trunc) {
        super(parent, trunc);
    }
    
    /**
     * Returns r dot grad(u), with any truncation applied.  Does not include
     * division by D, to avoid repeated multiplication of this term when summing
     * over all pairs.  Negation and division by D in most cases is required 
     * at some point when using this quantity.
     */
    public abstract double virial(AtomPair pair);
    
    public abstract double hyperVirial(AtomPair pair);
    
    public abstract Space.Vector gradient(AtomPair pair);
    
    /**
     * Integral used to evaluate correction to truncation of potential.
     */
    public abstract double integral(double rC);
    
}//end of Potential2Soft