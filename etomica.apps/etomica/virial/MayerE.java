package etomica.virial;

import etomica.AtomPair;
import etomica.potential.Potential2;
/**
 * @author kofke
 *
 * Simple e-bond, exp(-beta*u).
 */
public class MayerE extends MayerFunction {

	/**
	 * Constructor for MayerE.
	 */
	public MayerE(Potential2 potential) {
		super();
		this.potential = potential;
	}

	/**
	 * @see etomica.virial.MayerFunction#f(etomica.AtomPair, double)
	 */
	public double f(AtomPair pair, double beta) {
		return Math.exp(-beta*potential.energy(pair));
	}

	private final Potential2 potential;

}
