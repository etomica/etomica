package etomica.virial.cluster;

import etomica.Atom;
import etomica.virial.CoordinatePairSet;
import etomica.virial.MayerHardSphere;

/**
 * @author kofke
 *
 * Defines a cluster as the difference between one cluster and another with the
 * same structure but using hard-sphere f-bonds.
 */
public class ExcessHS extends Cluster {

	/**
	 * Constructor for ExcessHS.
	 */
	public ExcessHS(Cluster cluster) {
		this(cluster, 1.0);
	}
	public ExcessHS(Cluster cluster, double sigmaHS) {
		super(cluster.bondGroup()[0].f, cluster);
		clusterHS = new Cluster(new MayerHardSphere(sigmaHS),cluster);
	}

	Cluster clusterHS;
	/**
	 * @see etomica.virial.Cluster#value(etomica.Atom, etomica.Atom, etomica.virial.PairSet, double)
	 */
	public double value(Atom atom1, Atom atom2, CoordinatePairSet pairs, double beta) {
		return super.value(atom1, atom2, pairs, beta)-clusterHS.value(atom1, atom2, pairs, beta);
	}

	/**
	 * @see etomica.virial.Cluster#value(etomica.Atom, etomica.virial.PairSet, double)
	 */
	public double value(Atom atom, CoordinatePairSet pairs, double beta) {
		return super.value(atom, pairs, beta)-clusterHS.value(atom, pairs, beta);
	}

	/**
	 * @see etomica.virial.Cluster#value(etomica.virial.PairSet, double)
	 */
	public double value(CoordinatePairSet pairs, double beta) {
		return super.value(pairs, beta)-clusterHS.value(pairs,beta);
	}

}
