package etomica.virial.cluster;

import etomica.virial.Cluster;
import etomica.virial.MayerFunction;

/**
 * @author kofke
 *
 * A cluster formed as a simple chain of bonds.
 */
public class Chain extends Cluster {

	/**
	 * Constructor for Chain.
	 * @param n number of points in ring
	 * @param weight weight associated with cluster
	 * @param MayerFunction bond joining each pair in ring
	 */
	public Chain(int n, double weight, MayerFunction f) {
		super(n, weight, new Cluster.BondGroup(f, Standard.chain(n)));
	}
}
