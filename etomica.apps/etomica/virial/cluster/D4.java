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
		this(f, false);
	}
	public D4(MayerFunction f, boolean usePermutations) {
		super(4, new BondGroup[] {new BondGroup(f, Standard.D4)}, usePermutations);
	}
}
