package etomica.virial.cluster;

import etomica.virial.Cluster;

/**
 * @author kofke
 *
 * The virial cluster of 4 points joined as a simple ring.
 */
public final class D4 extends Cluster {
	public D4() {
		super(-3./8., new int[][] {{0,1},{0,3},{1,2},{2,3}});
	}
}
