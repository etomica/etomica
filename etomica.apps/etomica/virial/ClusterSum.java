package etomica.virial;

import etomica.Atom;

/**
 * @author kofke
 *
 * To change this generated comment edit the template variable "typecomment":
 * Window>Preferences>Java>Templates.
 * To enable and disable the creation of type comments go to
 * Window>Preferences>Java>Code Generation.
 */
public class ClusterSum implements ClusterAbstract {

	/**
	 * Constructor for ClusterSum.
	 */
	public ClusterSum(double weight, ClusterAbstract[] clusters) {
		super();
		this.weight = weight;
		nCluster = clusters.length;
		this.clusters = new ClusterAbstract[nCluster];
		nPoint = clusters[0].pointCount();
		for(int i=0; i<nCluster; i++) {
			this.clusters[i] = clusters[i];
			if(clusters[i].pointCount() != nPoint) throw new IllegalArgumentException("Attempt to construct ClusterSum with clusters having differing numbers of points");
		}
	}

	/**
	 * @see etomica.virial.ClusterValuable#pointCount()
	 */
	public int pointCount() {
		return nPoint;
	}

	/**
	 * @see etomica.virial.ClusterValuable#value(etomica.virial.PairSet, double)
	 */
	public double value(CoordinatePairSet cPairs, double beta) {
		double sum = 0.0;
		for(int i=0; i<nCluster; i++) {
			sum += clusters[i].weight() * clusters[i].value(cPairs, beta);
		}
		return sum;
	}

	/**
	 * @see etomica.virial.ClusterValuable#weight()
	 */
	public double weight() {
		return weight;
	}
	
	public ClusterAbstract[] cluster() {return clusters;}

	private final ClusterAbstract[] clusters;
	private final int nCluster, nPoint;
	final double weight;
}
