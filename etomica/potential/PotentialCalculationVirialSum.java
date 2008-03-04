package etomica.potential;

import etomica.api.IAtomSet;
import etomica.api.IPotential;
import etomica.atom.iterator.AtomsetIterator;

/**
 * Evaluates the virial summed over all iterated atoms.
 *
 * @author David Kofke
 */
public class PotentialCalculationVirialSum extends PotentialCalculation {
		
    /**
	 * Adds to the virial sum the energy values obtained from application of the given potential to the
	 * atoms produced by the given iterator.  Iterator is reset by method before beginning calculation.
	 */
	protected void doCalculation(AtomsetIterator iterator, IPotential potential) {
        if (!(potential instanceof PotentialSoft)) {
            return;
        }
		iterator.reset();
        for (IAtomSet atoms = iterator.next(); atoms !=null; atoms = iterator.next()) {
            sum += ((PotentialSoft)potential).virial(atoms);
		}
	}
	
	/**
	 * Sets the virial sum to zero, typically to begin a new virial-sum calculation.
	 * @return this instance, so the method can be called in-line as the instance is
	 * passed to the PotentialMaster.
	 */
	public PotentialCalculationVirialSum zeroSum() {
		sum = 0.0;
		return this;
	}

    private static final long serialVersionUID = 1L;

	/**
	 * Returns the current value of the energy sum.
	 */
	public double getSum() {return sum;}
	
	private double sum = 0.0;

 }//end VirialSum
