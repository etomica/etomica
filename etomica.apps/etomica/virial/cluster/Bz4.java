package etomica.virial.cluster;

import etomica.virial.Cluster;
import etomica.virial.ClusterAbstract;
import etomica.virial.ClusterSum;
import etomica.virial.MayerE;

/**
 * @author kofke
 *
 * Cluster sum for the 4th-order terms in the z (activity) expansion of the
 * pressure.
 */
public class Bz4 extends ClusterSum {

	/**
	 * Constructor for Bz4.
	 * @param weight
	 * @param clusters
	 */
	public Bz4(MayerE e) {
		super(new ClusterAbstract[] {
						new Cluster(4, new Cluster.BondGroup(e, Standard.D6)),
						new Cluster(4, new Cluster.BondGroup[] {new Cluster.BondGroup(e, Standard.C3)}, true),
						new Cluster(4, new Cluster.BondGroup[] {new Cluster.BondGroup(e, new int[][] {{0,1},{2,3}})}, true),
						new Cluster(4, new Cluster.BondGroup[] {new Cluster.BondGroup(e, new int[][] {{0,1}})}, true),
						new ClusterUnity(4)}, new double[] {+1.0, -4.0, -3.0, +12.0, -6.0} );
	}
}
