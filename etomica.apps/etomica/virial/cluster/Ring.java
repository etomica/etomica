package etomica.virial.cluster;

import etomica.virial.Cluster;
import etomica.virial.MayerFunction;

/**
 * @author kofke
 *
 * A cluster formed as a simple ring of equivalent bonds.
 */
public class Ring extends Cluster {

	/**
	 * Constructor for Ring.
	 * @param n number of points in ring
	 * @param weight weight associated with cluster
	 * @param MayerFunction bond joining each pair in ring
	 */
	public Ring(int n, double weight, MayerFunction f) {
		super(n, weight, new Cluster.BondGroup(f, Standard.ring(n)));
	}
}
