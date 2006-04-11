package etomica.potential;

import etomica.atom.AtomSet;
import etomica.atom.iterator.AtomsetIterator;
import etomica.util.Debug;

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
        if (iterator.hasNext() && iterator.peek() == null) {
            throw new RuntimeException("oops");
        }
		while(iterator.hasNext()) {
            AtomSet atoms = iterator.next();
            double e = potential.energy(atoms);
            if ((Debug.DEBUG_NOW || true) && e == Double.POSITIVE_INFINITY) {
                if (Debug.DEBUG_NOW) {
                    System.out.println("overlap for "+atoms);
                    potential.energy(atoms);
                }
                overlaps++;
            }
			sum += e;
		}
	}
	
	/**
	 * Sets the energy sum to zero, typically to begin a new energy-sum calculation.
	 * @return this instance, so the method can be called in-line as the instance is
	 * passed to the PotentialMaster.
	 */
	public PotentialCalculationEnergySum zeroSum() {
		sum = 0.0;
        overlaps = 0;
		return this;
	}

	/**
	 * Returns the current value of the energy sum.
	 */
	public double getSum() {
        if (overlaps > 0) {
            System.out.println(overlaps + " overlaps");
        }
        return sum;
    }
	
    private int overlaps = 0;
	private  double sum = 0.0;
        
}
