package etomica.virial;

import etomica.potential.Potential2Spherical;
import etomica.space.CoordinatePair;
/**
 * @author kofke
 *
 * Simple e-bond, exp(-beta*u).
 */
public class MayerE implements MayerFunction {

	/**
	 * Constructor for MayerE.
	 */
	public MayerE(Potential2Spherical potential) {
		this.potential = potential;
	}

	/**
	 * @see etomica.virial.MayerFunction#f(etomica.AtomPair, double)
	 */
	public double f(CoordinatePair cPair, double beta) {
		return Math.exp(-beta*potential.u(cPair.r2()));
	}

	private final Potential2Spherical potential;

}
