package etomica;

public interface Potential3Calculation extends PotentialCalculation {
    
    public void calculate(Atom3Iterator iterator, Potential3 potential);

    /*  contemplated for addition
    public Atom3Action getAtom3Calculation(Potential3 potential);
    */
}//end of Potential3Calculation
