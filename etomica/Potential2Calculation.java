package etomica;

public interface Potential2Calculation extends PotentialCalculation {
    
    public void calculate(AtomPairIterator iterator, Potential2 potential);

    /*  contemplated for addition
    public AtomPairAction getAtomPairCalculation(Potential2 potential);
    */
}//end of Potential2Calculation
