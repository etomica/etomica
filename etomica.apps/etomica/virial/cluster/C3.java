package etomica.virial.cluster;

import etomica.virial.Cluster;

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
	public C3() {
		super(-1./3.,new int[][] {{0,1},{0,2},{1,2}});
	}
}
