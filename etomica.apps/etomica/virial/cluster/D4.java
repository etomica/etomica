package etomica.virial.cluster;

import etomica.virial.Cluster;
import etomica.virial.MayerFunction;

/**
 * @author kofke
 *
 * The virial cluster of 4 points joined as a simple ring.
 */
public final class D4 extends Cluster {
	public D4(MayerFunction f) {
		super(4, -3./8., new BondGroup(f, Standard.D4));
	}
}
