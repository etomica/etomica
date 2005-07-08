package etomica.potential;

import etomica.AtomsetIterator;
import etomica.IteratorDirective;
import etomica.Potential;
import etomica.PotentialGroup;
import etomica.PotentialMaster;

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

public abstract class PotentialCalculation implements java.io.Serializable {
 	
	/**
	 * Method called to perform the calculation defined by this class.  If the given
	 * potential is an instance of PotentialGroup, the doCalculation action is to
	 * call the potential's calculate method, passing it the iterator, directive, and this
	 * instance. If the potential is not a PotentialGroup, then this class's
	 * doCalculation(AtomsetIterator, Potential) method is called with the given iterator
	 * and potential.
	 * 
	 * @param iterator
	 * @param id
	 * @param potential
	 */
	public void doCalculation(AtomsetIterator iterator, IteratorDirective id, Potential potential) {	
		if(potential instanceof PotentialGroup) {
			((PotentialGroup)potential).calculate(iterator, id, this);
		} else {
			doCalculation(iterator, potential);
		}
	}
	
	/**
	 * Method giving the specific calculation performed by this class.  Concrete subclasses
	 * are expected to invoke the iterator's reset() method before beginning iteration.
	 * @param iterator Iterator that has been conditioned to return the desired atom-sets
	 * for which the calculation is performed.
	 * @param potential The potential used to apply the action defined by this class.
	 */
	protected abstract void doCalculation(AtomsetIterator iterator, Potential potential);
	
}//end of PotentialCalculation
