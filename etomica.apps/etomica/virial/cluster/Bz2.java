package etomica.virial.cluster;

import etomica.virial.Cluster;
import etomica.virial.MayerFunction;

/**
 * @author kofke
 *
 * Cluster for the 2nd order terms in the z (activity) expansion of the
 * pressure.
 */
public class Bz2 extends Cluster {

	/**
	 * Constructor for Bz2.
	 * @param n
	 * @param weight
	 * @param bonds
	 */
	public Bz2(MayerFunction f) {
		super(2, new BondGroup(f, Standard.B2));
	}

}
