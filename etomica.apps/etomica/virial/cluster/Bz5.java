package etomica.virial.cluster;

import etomica.virial.Cluster;
import etomica.virial.ClusterSum;
import etomica.virial.MayerE;

/**
 * @author kofke
 *
 * Cluster sum for the 5th-order terms in the z (activity) expansion of the
 * pressure.
 */
public class Bz5 extends ClusterSum {

	/**
	 * Constructor for Bz5.
	 * @param clusters
	 */
	public Bz5(MayerE e) {
		super(new Cluster[] {
						new Cluster(5, new Cluster.BondGroup(e, Standard.full(5))),//Q5
						new Cluster(5, new Cluster.BondGroup[] {new Cluster.BondGroup(e, Standard.D6)}, true),//Q4 Q1
						new Cluster(5, new Cluster.BondGroup[] {new Cluster.BondGroup(e, new int[][] {{0,1},{0,2},{1,2},{3,4}})}, true),//Q3 Q2
						new Cluster(5, new Cluster.BondGroup[] {new Cluster.BondGroup(e, Standard.C3)}, true),//Q3 Q1^2
						new Cluster(5, new Cluster.BondGroup[] {new Cluster.BondGroup(e, new int[][] {{0,1},{2,3}})}, true),//Q2^2 Q1
						new Cluster(5, new Cluster.BondGroup[] {new Cluster.BondGroup(e, Standard.B2)}, true),//Q2 Q1^3
						new ClusterUnity(5)}, new double[] {+1.0, -5.0, -10.0, +20.0, +30.0, -60.0, +24.0} );//Q1^5
	}

}
