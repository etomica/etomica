package etomica.virial.overlap;

import etomica.MeterFunction;
import etomica.Simulation;
import etomica.units.Dimension;
import etomica.virial.Cluster;

/**
 * @author kofke
 *
 * To change this generated comment edit the template variable "typecomment":
 * Window>Preferences>Java>Templates.
 * To enable and disable the creation of type comments go to
 * Window>Preferences>Java>Code Generation.
 */
public class MeterOverlap extends MeterFunction implements etomica.DatumSource {

	/**
	 * Constructor for MeterOverlap.
	 * @param sim
	 */
	public MeterOverlap(Simulation sim, Cluster cluster0, Cluster cluster1) {
		super(sim);
		setX(-5, 5, 100);
		this.cluster0 = cluster0;
		this.cluster1 = cluster1;
		System.out.println("MeterOverlap constructor: "+x[nPoints/2]);
	}

	/**
	 * @see etomica.DataSource.X#getXDimension()
	 */
	public Dimension getXDimension() {
		return Dimension.NULL;
	}

	/**
	 * @see etomica.MeterFunction#currentValue()
	 */
	public double[] currentValue() {
		double v0 = Math.abs(cluster0.value(beta));
		double v1 = cluster1.value(beta);
		if(signPositive != (v1>0)) v1 = 0.0;
		else v1 = Math.abs(v1);
//		System.out.println(v1);
		for(int i=0; i<nPoints; i++) {
			y[i] = 1.0/v0/(1.0/v0 + Math.exp(x[i])/v1);
			y[i] = 1.0/(1.0 + Math.exp(x[i])*v0/v1);
		}
		return y;
	}
	
	public double value(etomica.DataSource.ValueType dummy) {
		return average()[nPoints/2];
	}

	/**
	 * @see etomica.MeterAbstract#getDimension()
	 */
	public Dimension getDimension() {
		return Dimension.NULL;
	}

	private Cluster cluster0, cluster1;
	double temperature, beta;
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

	private boolean signPositive;
	/**
	 * Returns the signPositive.
	 * @return boolean
	 */
	public boolean isSignPositive() {
		return signPositive;
	}

	/**
	 * Sets the signPositive.
	 * @param signPositive The signPositive to set
	 */
	public void setSignPositive(boolean signPositive) {
		this.signPositive = signPositive;
	}

}
