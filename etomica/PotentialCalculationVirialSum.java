package etomica;

/**
 * Evaluates the virial summed over all iterated atoms.
 *
 * @author David Kofke
 */
 
 /* History
  * 10/12/02 (DAK) new
  */
public class PotentialCalculationVirialSum implements PotentialCalculation.Sum, 
                                                        Potential2Calculation {
    protected double sum = 0.0;
        
    public PotentialCalculationVirialSum() {}
        
    public PotentialCalculation.Sum reset() {sum = 0.0; return this;}
    public double sum() {return sum;}
    
    //pair
    public void calculate(AtomPairIterator iterator, Potential2 potential) {
        if(!(potential instanceof Potential2Soft)) throw new RuntimeException("Error: PotentialCalculationVirialSum being used with potential that is not soft 2-body type");
        while(iterator.hasNext()) {
            sum += ((Potential2Soft)potential).virial(iterator.next());
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
            sum += ((Potential2Soft)potential).virial(pair);
        }
    }
        
}//end VirialSum
