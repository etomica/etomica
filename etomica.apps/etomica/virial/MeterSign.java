package etomica.virial;

import etomica.MeterScalar;
import etomica.SimulationElement;
import etomica.units.Dimension;

/**
 * @author kofke
 */

/* History
 * 08/21/03 (DAK) invoke resetPairs for pairSet in currentValue
 */
public class MeterSign extends MeterScalar {

	private Cluster cluster;
	private double temperature, beta;
	
	/**
	 * Constructor for MeterSign.
	 * @param parent
	 */
	public MeterSign(SimulationElement parent, Cluster cluster) {
		super(parent);
		setCluster(cluster);
		setTemperature(1.0); //almost always, temperature is irrelevant.  Sign of f does not normally depend on temperature.
	}

	/**
	 * @see etomica.MeterScalar#getData()
	 */
	public double getData() {
		return (cluster.value(((PhaseCluster)phase).getPairSet().resetPairs(), beta)>0) ? +1.0 : -1.0;
	}

	/**
	 * @see etomica.MeterAbstract#getDimension()
	 */
	public Dimension getDimension() {
		return null;
	}

	/**
	 * Returns the cluster.
	 * @return Cluster
	 */
	public Cluster getCluster() {
		return cluster;
	}

	/**
	 * Returns the temperature.
	 * @return double
	 */
	public double getTemperature() {
		return temperature;
	}

	/**
	 * Sets the cluster.
	 * @param cluster The cluster to set
	 */
	public void setCluster(Cluster cluster) {
		this.cluster = cluster;
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
