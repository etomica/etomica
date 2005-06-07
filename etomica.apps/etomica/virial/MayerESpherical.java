package etomica.virial;

import etomica.Space;
import etomica.potential.Potential2Spherical;
import etomica.space.CoordinatePair;
/**
 * @author kofke
 *
 * Simple e-bond, exp(-beta*u).
 */
public class MayerESpherical extends MayerFunctionSpherical {

	/**
	 * Constructor for MayerESpherical.
	 */
	public MayerESpherical(Space space, Potential2Spherical potential) {
        super(space);
		this.potential = potential;
	}

	/**
	 * @see etomica.virial.MayerFunctionSpherical#f(etomica.AtomPair, double)
	 */
	public double f(CoordinatePair cPair, double beta) {
		return Math.exp(-beta*potential.u(cPair.r2()));
	}

	private final Potential2Spherical potential;

}
