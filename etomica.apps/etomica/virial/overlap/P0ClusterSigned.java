package etomica.virial.overlap;

import etomica.IteratorDirective;
import etomica.Phase;
import etomica.PotentialCalculation;
import etomica.PotentialCalculationEnergySum;
import etomica.SimulationElement;
import etomica.virial.*;

/**
 * @author kofke
 *
 * Pair potential based on the value of a cluster, but giving an infinite energy
 * if the sign of the value is not positive or negative (depending on
 * specification given via setSign method.
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
	 * @see etomica.Potential0#energy(etomica.Phase)
	 */
	public void calculate(Phase phase, IteratorDirective id, PotentialCalculation pc) {
	   if(!enabled) return;
	   double beta = 1.0/phase.integrator().temperature();
	   PairSet pairs = ((PhaseCluster)phase).getPairSet();
	   g = cluster.value(pairs, beta);
//	   switch(id.atomCount()) {
//			case 0: g = cluster.value(pairs, beta); break;
//			case 1: g = cluster.value(id.atom1(), pairs, beta); break;
//			default: throw new RuntimeException();
//	   }
	   ((PotentialCalculationEnergySum)pc).set(this).actionPerformed(phase);
	}//end of calculate

	public double energy(Phase phase) {
		if(g == 0) return Double.POSITIVE_INFINITY;
//		double absg = (g>0) ? g : -g;
		if( signPositive != (g>0)) return Double.POSITIVE_INFINITY;
		else return -phase.integrator().temperature()*Math.log( (g>0) ? g : -g); //arg to log is abs(g)
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
