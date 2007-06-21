package etomica.potential;

import etomica.atom.AtomSet;
import etomica.atom.iterator.AtomsetIterator;

/**
 * Sums the number of iterates given by an iterator.  Useful as a debugging
 * tool.
 *
 * @author David Kofke
 */
public final class PotentialCalculationIterateSum extends PotentialCalculation {

    private static final long serialVersionUID = 1L;
    double sum = 0.0;
        
    public PotentialCalculationIterateSum() {}
        
    public PotentialCalculationIterateSum reset() {
    	sum = 0.0; 
    	return this;
    }
    
    public double getSum() {return sum;}
    
	public void doCalculation(AtomsetIterator iterator, IPotential potential) {
		iterator.reset();
        for (AtomSet atoms = iterator.next(); atoms != null; atoms = iterator.next()) sum++;
	}      
	
}//end PairSum
