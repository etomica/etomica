package etomica.virial;

import etomica.MeterScalar;
import etomica.SimulationElement;
import etomica.units.Dimension;

/**
 * @author kofke
 *
 * To change this generated comment edit the template variable "typecomment":
 * Window>Preferences>Java>Templates.
 * To enable and disable the creation of type comments go to
 * Window>Preferences>Java>Code Generation.
 */
public class MeterSign extends MeterScalar {

	private Cluster cluster;
	private double temperature, beta;
	private PairSet pairSet;
	
	/**
	 * Constructor for MeterSign.
	 * @param parent
	 */
	public MeterSign(SimulationElement parent, PairSet pairSet, Cluster cluster) {
		super(parent);
		setCluster(cluster);
		setPairSet(pairSet);
		setTemperature(1.0); //almost always, temperature is irrelevant.  Sign of f does not normally depend on temperature.
	}

	/**
	 * @see etomica.MeterScalar#currentValue()
	 */
	public double currentValue() {
		return (cluster.value(pairSet, beta)>0) ? +1.0 : -1.0;
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

	/**
	 * Returns the pairSet.
	 * @return PairSet
	 */
	public PairSet getPairSet() {
		return pairSet;
	}

	/**
	 * Sets the pairSet.
	 * @param pairSet The pairSet to set
	 */
	public void setPairSet(PairSet pairSet) {
		this.pairSet = pairSet;
	}

}
