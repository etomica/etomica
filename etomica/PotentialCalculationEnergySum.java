package etomica;

/**
 * Evaluates the energy summed over all iterated atoms.
 */
public class PotentialCalculationEnergySum implements PotentialCalculation.Sum, 
                                                        Potential0Calculation,
                                                        Potential1Calculation, 
                                                        Potential2Calculation {
    protected double sum = 0.0;
        
    public PotentialCalculationEnergySum() {}
        
    public PotentialCalculation.Sum reset() {sum = 0.0; return this;}
    public double sum() {return sum;}
    
    //zero-body
    public void calculate(Potential0 potential) {
        sum += potential.energy();
    }
        
    //atom
    public void calculate(AtomIterator iterator, Potential1 potential) {
        while(iterator.hasNext()) {
            sum += potential.energy(iterator.next());
            if(sum >= Double.MAX_VALUE) return;
        }//end while
    }//end of calculate

    //pair
    public void calculate(AtomPairIterator iterator, Potential2 potential) {
        while(iterator.hasNext()) {
            sum += potential.energy(iterator.next());
            if(sum >= Double.MAX_VALUE) return;
        }//end while
    }//end of calculate
    
    
    // the following methods/classes provide an alternate way to sum the energy
    // by providing actions to be passed to iterators, rather than passing
    // the iterators to this PotentialCalculation
    
    public AtomAction getAtomCalculation(Potential1 potential) {
        atomAction.potential = potential;
        return atomAction;
    }//end of getAtomCalculation
    
    public AtomPairAction getAtomPairCalculation(Potential2 potential) {
        atomPairAction.potential = potential;
        return atomPairAction;
    }//end of getAtomPairCalculation
    
    private final MyAtomAction atomAction = new MyAtomAction();
    private final class MyAtomAction extends AtomAction {
        Potential1 potential;
        public void actionPerformed(Atom atom) {
            sum += potential.energy(atom);
        }
    }
        
    private final MyAtomPairAction atomPairAction = new MyAtomPairAction();
    private final class MyAtomPairAction extends AtomPairAction {
        Potential2 potential;
        public void action(AtomPair pair) {
            sum += potential.energy(pair);
        }
    }
        
}//end EnergySum
