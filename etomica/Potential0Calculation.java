package etomica;

/**
 * Interface for a zero-body calculation, with a result that ususally
 * depends on more general features of a phase, such as its volume or
 * density.  Used for long-range correction contributions.
 *
 * @see Potential0Lrc
 * @author David Kofke
 */
public interface Potential0Calculation extends PotentialCalculation {
    
    public void calculate(Potential0 potential);
    
}
