package etomica;

public interface Potential1Soft {
    
    public double energy(Atom atom);
        
    public Space.Vector gradient(Atom atom);
    
}