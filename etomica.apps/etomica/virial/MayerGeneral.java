package etomica.virial;

import etomica.AtomPair;
import etomica.Potential2;

/**
 * @author kofke
 *
 * General Mayer function, which wraps the Mayer potential around an instance of
 * a Potential2 object.
 */
public class MayerGeneral extends MayerFunction {

	/**
	 * Constructor for MayerFunctionGeneral.
	 */
	public MayerGeneral(Potential2 potential) {
		super();
		this.potential = potential;
	}

	/**
	 * @see etomica.virial.MayerFunction#f(etomica.AtomPair, double)
	 */
	public double f(AtomPair pair, double beta) {
		return Math.exp(-beta*potential.energy(pair)) - 1.0;
//		if(pair.r2() > 1.0) return Math.exp(-beta*potential.energy(pair)) - 1.0;
//		else return Math.exp(-beta*potential.energy(pair));
	}

	private final Potential2 potential;
}
