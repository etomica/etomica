package etomica.potential;

import etomica.atom.AtomSet;

/**
 * Interface for a polarizable potential.  The interface includes only a method
 * that returns the polarization energy from the last call to energy(AtomSet).
 * 
 * @author Ken
 */
public interface PotentialPolarizable {

    /**
     * Returns the polarization energy from the last call to energy(AtomSet).
     * @return
     */
	public double getLastPolarizationEnergy();

	/**
	 * Returns the polarization contribution (non pairwise-additive) to the
	 * energy for the given AtomSet.
	 */
	public double getPolarizationEnergy(AtomSet atoms);
}
