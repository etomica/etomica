/*
 * History
 * Created on Aug 31, 2004 by kofke
 */
package etomica.potential;

import etomica.Atom;
import etomica.space.Vector;


/**
 * Methods for properties obtained for a soft, differentiable pair potential.
 *
 * @author David Kofke
 */
public interface PotentialSoft {
	   
	public double energy(Atom[] atoms);

	/**
	 * Returns the gradient of the potential, indicating how
	 * the energy would change as the position of the first atom
	 * is varied.
	 * @param atoms
	 * @return
	 */
	//TODO consider if implementation for pair describes energy change with movement of first or 2nd atom
	//TODO consider returning array of vectors, one for each atom
	public Vector gradient(Atom[] atoms);

}