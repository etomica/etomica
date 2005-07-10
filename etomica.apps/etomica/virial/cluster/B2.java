package etomica.virial.cluster;

import etomica.virial.MayerFunction;

/**
 * @author kofke
 *
 * The virial cluster of 2 points joined by a bond.
 */
public final class B2 extends Cluster {

	/**
	 * Constructor for B2.
	 * @param weight
	 * @param pairs
	 */
	public B2(MayerFunction f) {
		super(2, new BondGroup(f, Standard.B2));
	}
}
