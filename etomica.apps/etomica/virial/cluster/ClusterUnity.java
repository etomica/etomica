package etomica.virial.cluster;

import etomica.virial.Cluster;
import etomica.virial.CoordinatePairSet;

/**
 * @author kofke
 *
 * Cluster with value of unity regardless of configuration.
 */
public class ClusterUnity extends Cluster {

	public ClusterUnity(int nPoints) {
		super(1,new BondGroup[0]);
		this.nPoints = nPoints;
	}

	/**
	 * @see etomica.virial.ClusterValuable#pointCount()
	 */
	public int pointCount() {
		return nPoints;
	}

	/**
	 * @see etomica.virial.ClusterValuable#value(etomica.virial.PairSet, double)
	 */
	public double value(CoordinatePairSet pairs, double beta) {
		return 1.0;
	}

	/**
	 * @see etomica.virial.ClusterValuable#weight()
	 */
	public double weight() {
		return weight;
	}

	private double weight;
	int nPoints;
}
