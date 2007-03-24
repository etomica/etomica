package etomica.potential;

import etomica.atom.AtomSet;
import etomica.space.IVector;


/**
 * Methods for properties obtained for a soft, differentiable pair potential.
 *
 * @author David Kofke
 */
public interface PotentialSoft {

    public double virial(AtomSet atoms);

    /**
	 * Returns the gradient of the potential as it applies to each atom in the 
     * given AtomSet, indicating how the energy would change as the position of 
     * the first atom is varied.  The method is allowed to return an array of
     * Vectors with fewer elements than the number of atoms in the AtomSet.
	 * @param atoms
	 * @return
	 */
	public IVector[] gradient(AtomSet atoms);

}