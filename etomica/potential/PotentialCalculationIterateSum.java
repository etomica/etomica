package etomica.potential;

import etomica.atom.iterator.AtomsetIterator;

/**
 * Sums the number of iterates given by an iterator.  Useful as a debugging
 * tool.
 *
 * @author David Kofke
 */
public final class PotentialCalculationIterateSum extends PotentialCalculation {

	double sum = 0.0;
        
    public PotentialCalculationIterateSum() {}
        
    public PotentialCalculationIterateSum reset() {
    	sum = 0.0; 
    	return this;
    }
    
    public double getSum() {return sum;}
    
	public void doCalculation(AtomsetIterator iterator, Potential potential) {
		iterator.reset();
		while(iterator.hasNext()) sum++;
	}      
	
}//end PairSum
