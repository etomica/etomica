package etomica;

public interface Potential2Soft {
    
    public double energy(AtomPair pair);
    
    public double virial(AtomPair pair);
    
    public double hyperVirial(AtomPair pair);
    
    public Space.Vector gradient(AtomPair pair);
    
}