package etomica.virial.dos;

import etomica.*;
import etomica.virial.*;
import etomica.units.Dimension;

/**
 * @author kofke
 *
 * To change this generated comment edit the template variable "typecomment":
 * Window>Preferences>Java>Templates.
 * To enable and disable the creation of type comments go to
 * Window>Preferences>Java>Code Generation.
 */
public class MeterDOS extends MeterFunction {

	private P0Cluster p0;
	private Cluster cluster;
	private double dx;
	
	/**
	 * Constructor for MeterDOS.
	 * @param parent
	 * @param cluster
	 * @param temperature
	 */
	public MeterDOS(Simulation parent, Cluster cluster) {
		super(parent);
		setActive(true);
		setCluster(cluster);
		double xMax = Math.exp(3.0);
		setX(0.0, xMax, 100);
		dx = xMax/nPoints;
	}

	/**
	 * @see etomica.MeterScalar#currentValue()
	 */
	public double[] currentValue() {
		for(int i=0; i<nPoints; i++) y[i] = 0.0;
		int k = (int)(cluster.value(((PhaseCluster)phase).getPairSet(),1.0)/dx);
		if(k>=y.length) k = y.length-1;
		y[k] = 1.0;// /p0.pi(((PhaseCluster)phase));
		return y;
	}

	/**
	 * Returns the p0.
	 * @return P0Cluster
	 */
	public P0Cluster getP0() {
		return p0;
	}

	/**
	 * Sets the p0.
	 * @param p0 The p0 to set
	 */
	public void setP0(P0Cluster p0) {
		this.p0 = p0;
	}

	public Dimension getXDimension() {return Dimension.NULL;}
	
	public Dimension getDimension() {return Dimension.NULL;}

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

}
