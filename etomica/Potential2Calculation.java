package etomica;

public interface Potential2Calculation extends PotentialCalculation {
    
    public void calculate(AtomPairIterator iterator, Potential2 potential);
    
}
