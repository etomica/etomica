package etomica.potential;

import etomica.atom.iterator.AtomsetIterator;

/**
 * Evaluates the energy summed over all iterated atoms. Each call to doCalculate
 * accumulates the energy of the given potential applied to the atoms produced
 * by the given iterator. Each such call accumulates a single sum, which can
 * be accessed via the getSum() method.  Sum is re-zeroed only upon a call to the
 * zeroSum() method. 
 *
 * @author David Kofke
 */

/* History
 * 08/29/03 (DAK) added actionPerformed(AtomSet) method because method made
 * abstract in PotentialCalculation
 * 08/31/04 (DAK) overhauled with change in potentials/iterators
 */
public final class PotentialCalculationEnergySum extends PotentialCalculation {

	/**
	 * Adds to the energy sum the energy values obtained from application of the given potential to the
	 * atoms produced by the given iterator.  Iterator is reset by method before beginning calculation.
	 */
	protected void doCalculation(AtomsetIterator iterator, Potential potential) {
		iterator.reset();
		while(iterator.hasNext()) {
			sum += potential.energy(iterator.next());
		}
	}
	
	/**
	 * Sets the energy sum to zero, typically to begin a new energy-sum calculation.
	 * @return this instance, so the method can be called in-line as the instance is
	 * passed to the PotentialMaster.
	 */
	public PotentialCalculationEnergySum zeroSum() {
		sum = 0.0;
		return this;
	}

	/**
	 * Returns the current value of the energy sum.
	 */
	public double getSum() {return sum;}
	
	private  double sum = 0.0;
        
}//end EnergySum
