package etomica.virial;

import etomica.IteratorDirective;
import etomica.Phase;
import etomica.Potential0;
import etomica.PotentialCalculation;
import etomica.PotentialCalculationEnergySum;
import etomica.SimulationElement;

/**
 * @author David Kofke
 *
 * Pair potential given according to the Mayer bonds in a sum of cluster
 * integrals (using their absolute values).
 */
public class P0ClusterUmbrella extends P0Cluster {

	protected Cluster[] cluster;
	private int nCluster;
	/**
	 * Constructor for P0Cluster.
	 * @param parent
	 */
	public P0ClusterUmbrella(SimulationElement parent) {
		super(parent);
	}
	public P0ClusterUmbrella(SimulationElement parent, Cluster[] cluster) {
		this(parent);
		setCluster(cluster);
	}
	
	public void setCluster(Cluster cluster) {
		setCluster(new Cluster[] {cluster});
	}
	
	public void setCluster(Cluster[] cluster) {
		this.cluster = cluster;
		nCluster = cluster.length;
	}
	
	public double pi(PhaseCluster phase) {
		double beta = 1.0/phase.integrator().temperature();
		PairSet pairs = ((PhaseCluster)phase).getPairSet();
		double pi = 0.0;
		for(int i=0; i<nCluster; i++) {
			double value = cluster[i].value(pairs, beta);
			pi += (value>=0) ? value : -value;
		}
		return pi;
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
//			case 0: g = cluster.value(pairs, beta); break;
//			case 1: g = cluster.value(id.atom1(), pairs, beta); break;
//			default: throw new RuntimeException();
//	   }
//	   ((PotentialCalculationEnergySum)pc).set(this).actionPerformed(phase);
//	}//end of calculate

//	public double energy(Phase phase) {
//		double g = pi((PhaseCluster)phase);
//		return (g==1.0) ? 0.0 : -phase.integrator().temperature()*Math.log((g>0)?g:-g);
//	}

}
