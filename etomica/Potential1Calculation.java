package etomica;

public interface Potential1Calculation extends PotentialCalculation {
    
    public void calculate(AtomIterator iterator, Potential1 potential);
    
}