package etomica.virial.dos;

import etomica.virial.Cluster;
import etomica.virial.PairSet;

/**
 * @author kofke
 *
 * Evaluates to zero if one of the molecules is not within a cutoff distance of
 * one of the other ones.
 */
public class ClusterGrouped extends Cluster {

	private double rCutoff, rCutoff2;
	/**
	 * Constructor for ClusterGrouped.
	 * @param n
	 * @param weight
	 * @param bonds
	 */
	public ClusterGrouped(double rCutoff) {
		super(0, 1.0, BondGroup.NULL);
		setRCutoff(rCutoff);
	}

	public double value(PairSet pairs, double beta) {
		double p = 1.0;
		pairs.resetPairs();
		boolean in01 = pairs.getPair(0,1).r2() < rCutoff2;
		boolean in02 = pairs.getPair(0,2).r2() < rCutoff2;
		boolean in12 = pairs.getPair(1,2).r2() < rCutoff2;		
		if( (in01 || in02) && (in01 || in12) && (in02 || in12)) return 1.0;
		else return 0.0;
	}

	/**
	 * Returns the rCutoff.
	 * @return double
	 */
	public double getRCutoff() {
		return rCutoff;
	}

	/**
	 * Sets the rCutoff.
	 * @param rCutoff The rCutoff to set
	 */
	public void setRCutoff(double rCutoff) {
		this.rCutoff = rCutoff;
		rCutoff2 = rCutoff * rCutoff;
	}

}
