package etomica;

public interface Potential1Calculation extends PotentialCalculation {
    
    public void calculate(AtomIterator iterator, Potential1 potential);

    /* contemplated for addition
    public AtomAction getAtomCalculation(Potential1 potential);
    */
    
}