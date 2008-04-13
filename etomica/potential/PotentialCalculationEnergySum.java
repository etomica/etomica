package etomica.potential;

import etomica.api.IAtomSet;
import etomica.api.IPotential;

/**
 * Evaluates the energy summed over all iterated atoms. Each call to doCalculate
 * accumulates the energy of the given potential applied to the atoms produced
 * by the given iterator. Each such call accumulates a single sum, which can
 * be accessed via the getSum() method.  Sum is re-zeroed only upon a call to the
 * zeroSum() method. 
 *
 * @author David Kofke
 */
public final class PotentialCalculationEnergySum implements PotentialCalculation, java.io.Serializable {

    /**
	 * Adds to the energy sum the energy values obtained from application of the given potential to the
	 * atoms produced by the given iterator.  Iterator is reset by method before beginning calculation.
	 */
	public void doCalculation(IAtomSet atoms, IPotential potential) {
	    sum += potential.energy(atoms);
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
	public double getSum() {
        return sum;
    }
	
    private static final long serialVersionUID = 1L;
	private  double sum = 0.0;
}
