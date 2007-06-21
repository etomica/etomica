package etomica.potential;

/**
 * interface for spherical 2-body potentials
 */
public interface Potential2Spherical extends IPotential {
	/**
	 * The pair energy u(r^2) with no truncation applied.
	 * @param the square of the distance between the particles.
	 */
	public double u(double r2);

}