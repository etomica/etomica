package etomica.virial;

/**
 * @author kofke
 *
 * To change this generated comment edit the template variable "typecomment":
 * Window>Preferences>Java>Templates.
 * To enable and disable the creation of type comments go to
 * Window>Preferences>Java>Code Generation.
 */
public interface ClusterAbstract {

	/**
	 * Number of points in the cluster.
	 * @return int
	 */
	public int pointCount();
	
	/**
	 * Value of this cluster for the given pairset at the specified reciprocal
	 * temperature.
	 * @param pairs
	 * @param beta Reciprocal temperature
	 * @return double
	 */
	public double value(PairSet pairs, double beta);

	/**
	 * Weight associated with this cluster.
	 * @return double
	 */
	public double weight();
}
