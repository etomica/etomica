package etomica.virial.overlap;

import etomica.MeterFunctionGroup;
import etomica.Simulation;
import etomica.units.Dimension;
import etomica.virial.Cluster;
import etomica.virial.PairSet;
import etomica.virial.PhaseCluster;

/**
 * @author kofke
 *
 * Overlap-sampling evaluation of ratio of clusters, simulating the reference
 * cluster.
 */

/* History
 * 08/21/03 (DAK) invoke resetPairs on pairSet in currentValue method
 */
 
public class MeterOverlapReference extends MeterFunctionGroup {

	/**
	 * Constructor for MeterOverlap.
	 * @param sim
	 */
	public MeterOverlapReference(Simulation sim, Cluster cluster0, Cluster cluster1) {
		super(sim, 2);
		setX(-5, 5, 1000);
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
	 * @see etomica.MeterFunction#getData()
	 */
	public void updateValues() {
		PairSet pairSet = ((PhaseCluster)phase).getPairSet().resetPairs();
		double v0 = Math.abs(cluster0.value(pairSet,beta));
		double v1 = cluster1.value(pairSet,beta);
		boolean positive = (v1>0);
		if(!positive) v1 *= -1;//abs(v1)
		int i0 = (positive) ? 0 : 1;//index for average adding zero
		int i1 = 1-i0; //index for average adding value; 1 or 0 for i0 = 0 or 1
		for(int i=0; i<nPoints; i++) {
			y[i0][i] = 0.0;
			y[i1][i] = 1.0/(1.0 + expX[i]*v0/v1);
		}
	}

	/**
	 * @see etomica.MeterFunction#setX(double, double, int)
	 */
	public void setX(double min, double max, int n) {
		super.setX(min, max, n);
		expX = new double[n];
		for(int i=0; i<n; i++) {
			expX[i] = Math.exp(meters[0].X()[i]);
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
	double[] expX;
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
