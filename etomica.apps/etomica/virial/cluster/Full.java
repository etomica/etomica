package etomica.virial.cluster;

import etomica.virial.Cluster;
import etomica.virial.MayerFunction;

/**
 * @author kofke
 *
 * A cluster formed with a bond between every pair.
 */
public class Full extends Cluster {

	/**
	 * Constructor for a fully joined cluster
	 * @param n number of points in ring
	 * @param weight weight associated with cluster
	 * @param MayerFunction bond joining each pair in ring
	 */
	public Full(int n, double weight, MayerFunction f) {
		super(n, weight, new Cluster.BondGroup(f, Standard.full(n)));
	}
}
