package etomica.virial;

import etomica.Phase;
import etomica.SimulationElement;
import etomica.potential.Potential0;

/**
 * @author David Kofke
 *
 * Pair potential given according to the Mayer bonds in a cluster integral.
 * Does not require that the value of the cluster is non-negative.
 */

/* History
 * 08/20/03 (DAK) small changes to energy method (check for g = 0; abs(g)->g in
 * log argument
 * 08/21/03 (DAK) invoke resetPairs for pairSet in pi method
 * 12/16/03 (DAK) added field to hold indication of whether cluster giving value
 * of potential is positive or negative
 */
public class P0Cluster extends Potential0 {

	protected ClusterAbstract cluster;
	private boolean isPositive;
//	protected double g = 0;
	/**
	 * Constructor for P0Cluster.
	 * @param parent
	 */
	public P0Cluster(SimulationElement parent) {
		super(parent);
	}
	public P0Cluster(SimulationElement parent, ClusterAbstract cluster) {
		this(parent);
		setCluster(cluster);
	}
	
	public void setCluster(ClusterAbstract cluster) {
		this.cluster = cluster;
	}
	/**
	 * Returns the cluster.
	 * @return ClusterAbstract
	 */
	public ClusterAbstract getCluster() {
		return cluster;
	}
	
	public double pi(PhaseCluster phase) {
		double pi = cluster.value(phase.getPairSet().resetPairs(), 1.0/phase.integrator().temperature());
		isPositive = (pi > 0);
		return isPositive ? pi : -pi;
	}
//	public double pi(PairSet pairSet, double beta) {
//		double pi = cluster.value(pairSet, beta);
//		return (pi>0) ? pi : -pi;
//	}

//	public void calculate(Phase phase, IteratorDirective id, PotentialCalculation pc) {
//	   if(!enabled) return;
//	   double beta = 1.0/phase.integrator().temperature();
//	   PairSet pairs = ((PhaseCluster)phase).getPairSet();
//	   switch(id.atomCount()) {
//	   		case 0: g = cluster.value(pairs, beta); break;
//	   		case 1: g = cluster.value(id.atom1(), pairs, beta); break;
//	   		default: throw new RuntimeException();
//	   }
//	   ((PotentialCalculationEnergySum)pc).set(this).actionPerformed(phase);
//	}//end of calculate

	public double energy(Phase phase) {
		double g = pi((PhaseCluster)phase);
		return (g==1.0) ? 0.0 : ( (g==0) ? Double.POSITIVE_INFINITY : -phase.integrator().temperature()*Math.log(g) );
	}

	/**
	 * Returns true if the cluster value calculated during most recent call to
	 * the pi method was positive; returns false otherwise.
	 * @return boolean
	 */
	public boolean mostRecentIsPositive() {return isPositive;}

}
