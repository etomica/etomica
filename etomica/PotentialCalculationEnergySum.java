package etomica;

/**
 * Evaluates the energy summed over all iterated atoms.
 */
public final class PotentialCalculationEnergySum implements PotentialCalculation.Sum, Potential1Calculation, Potential2Calculation {
        
    private double sum = 0.0;
        
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
}//end EnergySum
