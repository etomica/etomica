package etomica.virial;


public interface ClusterAbstract {

	/**
	 * Number of points in the cluster.
	 * @return int
	 */
	public int pointCount();
	
	/**
	 * Value of this cluster for the given pairset at the specified reciprocal
	 * temperature.
	 */
	public double value(CoordinatePairSet cPairs, double beta);

}
