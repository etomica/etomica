package etomica.virial.cluster;

import etomica.virial.Cluster;

/**
 * @author kofke
 *
 * The virial cluster of 4 points joined as a simple ring, with one bond
 * joining an opposite pair.
 */
public final class D5 extends Cluster {
	public D5() {
		super(-3./4., new int[][] {{0,1},{0,2},{0,3},{1,2},{2,3}});
	}
}
