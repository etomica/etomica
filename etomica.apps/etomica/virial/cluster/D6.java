package etomica.virial.cluster;

import etomica.virial.Cluster;
import etomica.virial.MayerFunction;

/**
 * @author kofke
 *
 * The virial cluster of 4 points with all pairs joined.
 */
public final class D6 extends Cluster {
	public D6(MayerFunction f) {
		super(4, -1./8., new BondGroup(f, Standard.D6));
	}
}
