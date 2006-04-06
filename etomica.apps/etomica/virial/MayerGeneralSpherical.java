package etomica.virial;

import etomica.potential.Potential2Spherical;
import etomica.space.Space;

/**
 * @author kofke
 *
 * General Mayer function, which wraps the Mayer potential around an instance of
 * a Potential2 object.
 */
public class MayerGeneralSpherical extends MayerFunctionSpherical {

	/**
	 * Constructor Mayer function using given potential.
	 */
	public MayerGeneralSpherical(Space space, Potential2Spherical potential) {
        super(space);
		this.potential = potential;
	}

	public double f(double r2, double beta) {
		return Math.exp(-beta*potential.u(r2)) - 1.0;
	}

	private final Potential2Spherical potential;
}
