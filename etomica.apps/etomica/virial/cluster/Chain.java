package etomica.virial.cluster;

import etomica.virial.ClusterBonds;

/**
 * @author kofke
 *
 * A cluster formed as a simple chain of bonds.
 */
public class Chain extends ClusterBonds {

	/**
	 * Constructor for Chain.
	 * @param n number of points in ring
	 * @param weight weight associated with cluster
	 * @param MayerFunction bond joining each pair in ring
	 */
	public Chain(int n) {
		super(n, new int[][][] {Standard.chain(n)});
	}
}
