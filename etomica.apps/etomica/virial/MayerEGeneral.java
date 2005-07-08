package etomica.virial;

import etomica.AtomSet;
import etomica.Potential;

/**
 * @author kofke
 *
 * General Mayer function, which wraps the Mayer potential around an instance of
 * a Potential2 object.
 */
public class MayerEGeneral implements MayerFunction, java.io.Serializable {

	/**
	 * Constructor Mayer function using given potential.
	 */
	public MayerEGeneral(Potential potential) {
		this.potential = potential;
	}

	public double f(AtomSet pair, double beta) {
		return Math.exp(-beta*potential.energy(pair));
	}

	private final Potential potential;
}
