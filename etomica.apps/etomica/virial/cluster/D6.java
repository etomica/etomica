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
		this(f, false);
	}
	/**
	 * Ignores usePermutations argument, as this cluster has only one unique
	 * permutation.
	 */
	public D6(MayerFunction f, boolean usePermutations) {
		super(4, -1./8., new BondGroup[] {new BondGroup(f, Standard.D6)}, false);
	}
}
