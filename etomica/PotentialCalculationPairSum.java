package etomica;

/**
 * Sums the number of pairs returned by the iterator for the calculation.  
 * Useful as a debugging tool.
 *
 * @author David Kofke
 */
public class PotentialCalculationPairSum implements PotentialCalculation.Sum, 
                                                        Potential2Calculation {
    protected double sum = 0.0;
        
    public PotentialCalculationPairSum() {}
        
    public PotentialCalculation.Sum reset() {sum = 0.0; return this;}
    public double sum() {return sum;}
    
    //pair
    public void calculate(AtomPairIterator iterator, Potential2 potential) {
        while(iterator.hasNext()) {
            iterator.next();
            sum += 1;
        }//end while
    }//end of calculate
    
    // the following methods/classes provide an alternate way to sum the energy
    // by providing actions to be passed to iterators, rather than passing
    // the iterators to this PotentialCalculation
    
    public AtomPairAction getAtomPairCalculation(Potential2 potential) {
        atomPairAction.potential = potential;
        return atomPairAction;
    }//end of getAtomPairCalculation
    
    private final MyAtomPairAction atomPairAction = new MyAtomPairAction();
    private final class MyAtomPairAction extends AtomPairAction {
        Potential2 potential;
        public void action(AtomPair pair) {
            sum += 1;
        }
    }
        
}//end PairSum
