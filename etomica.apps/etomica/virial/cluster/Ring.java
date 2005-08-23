package etomica.virial.cluster;

import etomica.virial.ClusterBonds;

/**
 * @author kofke
 *
 * A cluster formed as a simple ring of equivalent bonds.
 */
public class Ring extends ClusterBonds {

	/**
	 * Constructor for Ring.
	 * @param n number of points in ring
	 * @param weight weight associated with cluster
	 * @param MayerFunction bond joining each pair in ring
	 */
	public Ring(int n) {
		super(n, new int[][][] {Standard.ring(n)});
	}
}
