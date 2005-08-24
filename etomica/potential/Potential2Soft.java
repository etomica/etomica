package etomica.potential;

import etomica.atom.AtomSet;

/**
 * Methods for properties obtained for a soft, differentiable pair potential.
 *
 * @author David Kofke
 */

/*
 * History
 * Created on Jan 21, 2005 by kofke
 */

public interface Potential2Soft extends PotentialSoft {
    
	public double hyperVirial(AtomSet pair);

	public double virial(AtomSet pair);
    
	/**
	 * Integral used to evaluate correction to truncation of potential.
	 */
	public abstract double integral(double rC);
}