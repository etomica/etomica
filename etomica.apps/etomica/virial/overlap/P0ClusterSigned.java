package etomica.virial.overlap;

import etomica.SimulationElement;
import etomica.virial.*;

/**
 * @author kofke
 *
 * Pair potential based on the value of a cluster, but giving an infinite energy
 * if the sign of the value is not positive or negative (depending on
 * specification given via setSign method.
 */

/* History
 * 08/20/03 (DAK) modified to make compatible with earlier revisions to
 * P0Cluster
 * 08/21/03 (DAK) invoke resetPairs for pairSet in pi method
 */
public class P0ClusterSigned extends P0Cluster {

	/**
	 * Constructor for P2ClusterSigned.
	 * @param parent
	 * @param pairs
	 */
	public P0ClusterSigned(SimulationElement parent) {
		super(parent);
		setSignPositive(true);
	}
	public P0ClusterSigned(SimulationElement parent, Cluster cluster) {
		super(parent, cluster);
		setSignPositive(true);
	}

	public boolean signPositive;
	
	/**
	 * Same as superclass method, but returns zero if sign of cluster is not in
	 * accord with signPositive value.  That is, returns zero if 
	 * [signPositive != (value > 0)]
	 * @see etomica.virial.P0Cluster#pi(PhaseCluster)
	 */
	public double pi(PhaseCluster phase) {
		double pi = cluster.value(phase.getPairSet().resetPairs(), 1.0/phase.integrator().temperature());
		if(signPositive != (pi>0)) return 0.0;
		else return (pi>0) ? pi : -pi;
	}
	
	/**
	 * Sets the sign.
	 * @param sign The sign to set
	 */
	public void setSignPositive(boolean signPositive) {
		this.signPositive = signPositive;
	}

	/**
	 * Returns the signPositive.
	 * @return boolean
	 */
	public boolean isSignPositive() {
		return signPositive;
	}

}
