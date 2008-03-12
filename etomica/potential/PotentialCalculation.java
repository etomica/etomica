package etomica.potential;

import etomica.api.IAtomSet;
import etomica.api.IPotential;

/**
 * Class defining a particular action to be performed on a set of atoms using an
 * arbitrary potential.  Examples of such actions are summing the energy, 
 * computing forces on atoms, determining collision times, etc.
 * Concrete subclasses define these actions through implementation of the 
 * doCalculation(AtomsetIterator, Potential) method, which should be set up
 * to loop through the iterates provided by the iterator (after resetting it)
 * and perform the defined calculation using the given potential.
 *
 * @see PotentialMaster
 * @see PotentialGroup
 */

public interface PotentialCalculation {
 	
	/**
	 * Method giving the specific calculation performed by this class.  Concrete subclasses
	 * are expected to invoke the iterator's reset() method before beginning iteration.
	 * @param iterator Iterator that has been conditioned to return the desired atom-sets
	 * for which the calculation is performed.
	 * @param potential The potential used to apply the action defined by this class.
	 */
	public void doCalculation(IAtomSet atoms, IPotential potential);
	
}
