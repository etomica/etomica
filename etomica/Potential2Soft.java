package etomica;

/**
 * Methods for properties obtained for a soft, differentiable pair potential.
 *
 * @author David Kofke
 */

public interface Potential2Soft {
    
    public double energy(AtomPair pair);
    
    public double virial(AtomPair pair);
    
    public double hyperVirial(AtomPair pair);
    
    public Space.Vector gradient(AtomPair pair);
    
    /**
     * Integral used to evaluate correction to truncation of potential.
     */
    public double integral(double rC);
    
}