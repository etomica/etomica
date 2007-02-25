package etomica.potential;

import etomica.atom.AtomSet;

/**
 * Methods for properties obtained for a soft, differentiable pair potential.
 *
 * @author David Kofke
 */
public interface Potential2Soft extends PotentialSoft {
    
	public double hyperVirial(AtomSet pair);

	/**
	 * Integral used to evaluate correction to truncation of potential.
	 */
	public abstract double integral(double rC);
}