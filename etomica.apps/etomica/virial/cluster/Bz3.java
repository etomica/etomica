package etomica.virial.cluster;

import etomica.virial.Cluster;
import etomica.virial.ClusterAbstract;
import etomica.virial.ClusterSum;
import etomica.virial.MayerE;

/**
 * @author kofke
 *
 * Cluster sum for the 3rd-order terms in the z (activity) expansion of the
 * pressure.
 */
public class Bz3 extends ClusterSum {

	/**
	 * Constructor for B3z.
	 * @param weight
	 * @param clusters
	 */
	public Bz3(MayerE e) {
		super(1./6., new ClusterAbstract[] {
						new Cluster(3, +1.0, new Cluster.BondGroup(e, Standard.C3)),
						new Cluster(3, -3.0, new Cluster.BondGroup[] {new Cluster.BondGroup(e, Standard.B2)}, true),
						new ClusterUnity(3, +2.0)} );
	}

}
