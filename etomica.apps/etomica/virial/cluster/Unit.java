package etomica.virial.cluster;

import etomica.virial.ClusterAbstract;
import etomica.virial.PairSet;

/**
 * @author kofke
 *
 * Cluster with value of unity regardless of configuration.
 */
public class Unit implements ClusterAbstract {

	/**
	 * Constructor for Unit, using default weight of 1.0
	 */
	public Unit(int nPoints) {
		this(nPoints, 1.0);
	}
	
	public Unit(int nPoints, double weight) {
		super();
		this.nPoints = nPoints;
		this.weight = weight;
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
	public double value(PairSet pairs, double beta) {
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
