package etomica.virial.overlap;

import etomica.MeterFunctionGroup;
import etomica.Simulation;
import etomica.units.Dimension;
import etomica.virial.Cluster;

/**
 * @author kofke
 *
 * Overlap-sampling evaluation of ratio of clusters, simulating the reference
 * cluster.
 */
public class MeterOverlapReference extends MeterFunctionGroup {

	/**
	 * Constructor for MeterOverlap.
	 * @param sim
	 */
	public MeterOverlapReference(Simulation sim, Cluster cluster0, Cluster cluster1) {
		super(sim, 2);
		setX(-5, 5, 100);
		this.cluster0 = cluster0;
		this.cluster1 = cluster1;
		allMeters()[0].setLabel("Negative");
		allMeters()[1].setLabel("Positive");
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
	public void updateValues() {
		double v0 = Math.abs(cluster0.value(beta));
		double v1 = cluster1.value(beta);
		boolean positive = (v1>0);
		if(!positive) v1 *= -1;//abs(v1)
		int i0 = (positive) ? 0 : 1;//index for average adding zero
		int i1 = 1-i0; //index for average adding value; 1 or 0 for i0 = 0 or 1
		for(int i=0; i<nPoints; i++) {
			y[i0][i] = 0.0;
			y[i1][i] = 1.0/(1.0 + Math.exp(meters[i1].X()[i])*v0/v1);
		}
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

}
