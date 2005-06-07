package etomica.virial;

import etomica.Phase;
import etomica.data.meter.MeterScalar;
import etomica.units.Dimension;

/**
 * @author kofke
 *
 * Returns the value of a cluster.
 */

/* History
 * 08/21/03 (DAK) invoke resetPairs for pairSet in currentValue
 */
public class MeterCluster extends MeterScalar {

	private Cluster cluster;
	private double temperature, beta;
	
	/**
	 * Constructor for MeterCluster.
	 * @param parent
	 */
	public MeterCluster(SimulationElement parent, Cluster cluster, double temperature) {
		super(parent);
		setCluster(cluster);
		setTemperature(temperature);
	}

	/**
	 * @see etomica.MeterScalar#getData()
	 */
	public double getDataAsScalar(Phase p) {
		return cluster.value(((PhaseCluster)p).getCPairSet(), beta);
	}

	/**
	 * @see etomica.MeterAbstract#getDimension()
	 */
	public Dimension getDimension() {
		return Dimension.NULL;
	}

	/**
	 * Returns the cluster.
	 * @return Cluster
	 */
	public Cluster getCluster() {
		return cluster;
	}

	/**
	 * Sets the cluster.
	 * @param cluster The cluster to set
	 */
	public void setCluster(Cluster cluster) {
		this.cluster = cluster;
	}

	/**
	 * Returns the temperature.
	 * @return double
	 */
	public double getTemperature() {
		return temperature;
	}

	/**
	 * Sets the temperature.
	 * @param temperature The temperature to set
	 */
	public void setTemperature(double temperature) {
		this.temperature = temperature;
		beta = 1.0/temperature;
	}
}
