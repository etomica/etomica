package etomica.virial;

import etomica.Atom;

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
	
    /**
     * Weight coefficient of the given cluster.
     * TODO: move this into ClusterSum, which is the only thing that uses it
     */
	public double weight();
	
}
