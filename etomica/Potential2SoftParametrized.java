package etomica;

public class Potential2SoftParametrized extends Potential2Soft {
    
    Potential2Soft potential;
    
    public Potential2SoftParametrized(PotentialGroup parent, Potential2Soft potential) {
        this.potential = potential;
    }
    
    public double energy(AtomPair pair) { }
        
    
    public double virial(AtomPair pair);
    
    public double hyperVirial(AtomPair pair);
    
    public Space.Vector gradient(AtomPair pair);
    
    /**
     * Integral used to evaluate correction to truncation of potential.
     */
    public double integral(double rC);
    
    
}