package etomica.virial.cluster;

import etomica.virial.Cluster;
import etomica.virial.MayerFunction;

/**
 * @author kofke
 *
 * The virial cluster of 3 points joined as a simple ring.
 */
public final class C3 extends Cluster {

	/**
	 * Constructor for C3.
	 * @param weight
	 * @param pairs
	 */
	public C3(MayerFunction f) {
		super(3, new BondGroup(f, Standard.C3));
	}
}
