package etomica;

/**
 * Evaluates the energy summed over all iterated atoms.
 */
public class PotentialCalculationEnergySum implements PotentialCalculation.Sum, Potential1Calculation, Potential2Calculation {
        
    protected double sum = 0.0;
        
    public PotentialCalculationEnergySum() {}
        
    public PotentialCalculation.Sum reset() {sum = 0.0; return this;}
    public double sum() {return sum;}
        
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
