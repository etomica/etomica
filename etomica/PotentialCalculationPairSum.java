package etomica;

/**
 * Sums the number of iterates given by an iterator.  Useful as a debugging
 * tool.
 *
 * @author David Kofke
 */
public class PotentialCalculationPairSum extends PotentialCalculation 
											implements PotentialCalculation.Summable {
    protected double sum = 0.0;
        
    public PotentialCalculationPairSum() {}
        
    public PotentialCalculation.Summable reset() {sum = 0.0; return this;}
    public double sum() {return sum;}
    
    public void actionPerformed(Atom atom) {
    	sum += 1;
    }
    public void actionPerformed(AtomPair pair) {
        sum += 1;
    }//end of calculate
    public void actionPerformed(Atom3 atom3) {
    	sum += 1;
    }    
        
}//end PairSum
