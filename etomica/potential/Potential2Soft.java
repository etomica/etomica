package etomica.potential;

import etomica.Atom;



/*
 * History
 * Created on Jan 21, 2005 by kofke
 */
/**
 * Methods for properties obtained for a soft, differentiable pair potential.
 *
 * @author David Kofke
 */
public interface Potential2Soft extends PotentialSoft {
	public double hyperVirial(Atom[] pair);

	public double virial(Atom[] pair);
    
	/**
	 * Integral used to evaluate correction to truncation of potential.
	 */
	public abstract double integral(double rC);
}