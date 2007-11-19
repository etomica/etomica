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
	public double integral(double rC);

    /**
     * The pair energy u(r^2).  Anisotropic potentials return an (Boltzmann
     * weighted) orientationally averaged energy.
     * @param the square of the distance between the particles.
     */
	public double u(double r2);
}