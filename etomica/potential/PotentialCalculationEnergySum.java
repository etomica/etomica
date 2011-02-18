package etomica.potential;

import etomica.api.IAtomList;
import etomica.api.IMoleculeList;
import etomica.api.IPotentialAtomic;
import etomica.api.IPotentialMolecular;

/**
 * Evaluates the energy summed over all iterated atoms. Each call to doCalculate
 * accumulates the energy of the given potential applied to the atoms produced
 * by the given iterator. Each such call accumulates a single sum, which can
 * be accessed via the getSum() method.  Sum is re-zeroed only upon a call to the
 * zeroSum() method. 
 *
 * @author David Kofke
 */
public class PotentialCalculationEnergySum implements PotentialCalculation, PotentialCalculationMolecular, java.io.Serializable {

    /**
	 * Adds to the energy sum the energy values obtained from application of the given potential to the
	 * atoms produced by the given iterator.  Iterator is reset by method before beginning calculation.
	 */
	public void doCalculation(IAtomList atoms, IPotentialAtomic potential) {
	    sum += potential.energy(atoms);
	}
	
    /**
     * Adds to the energy sum the energy values obtained from application of the given potential to the
     * atoms produced by the given iterator.  Iterator is reset by method before beginning calculation.
     */
    public void doCalculation(IMoleculeList atoms, IPotentialMolecular potential) {
        sum += potential.energy(atoms);
    }
    
	/**
	 * Sets the energy sum to zero, typically to begin a new energy-sum calculation.
	 */
	public void zeroSum() {
		sum = 0.0;
	}

	/**
	 * Returns the current value of the energy sum.
	 */
	public double getSum() {
        return sum;
    }
	
    private static final long serialVersionUID = 1L;
	protected double sum = 0.0;
}
