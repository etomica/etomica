package etomica.virial.cluster;

import etomica.virial.Cluster;

/**
 * @author kofke
 *
 * The virial cluster of 4 points with all pairs joined.
 */
public final class D6 extends Cluster {
	public D6() {
		super(-1./8., new int[][] {{0,1},{0,2},{0,3},{1,2},{1,3},{2,3}});
	}
}
