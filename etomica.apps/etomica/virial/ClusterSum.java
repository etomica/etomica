package etomica.virial;

public class ClusterSum implements ClusterAbstract {

	/**
	 * Constructor for ClusterSum.
	 */
	public ClusterSum(ClusterAbstract[] subClusters, double[] subClusterWeights) {
        if (subClusterWeights.length != subClusters.length) throw new IllegalArgumentException("number of clusters and weights must be the same");
		clusters = new ClusterAbstract[subClusters.length];
        clusterWeights = subClusterWeights;
		for(int i=0; i<clusters.length; i++) {
			clusters[i] = subClusters[i];
			if(clusters[i].pointCount() != clusters[0].pointCount()) throw new IllegalArgumentException("Attempt to construct ClusterSum with clusters having differing numbers of points");
		}
	}

	/**
	 * @see etomica.virial.ClusterValuable#pointCount()
	 */
	public int pointCount() {
		return clusters[0].pointCount();
	}

	/**
	 * @see etomica.virial.ClusterValuable#value(etomica.virial.PairSet, double)
	 */
	public double value(CoordinatePairSet cPairs, double beta) {
		double sum = 0.0;
		for(int i=0; i<clusters.length; i++) {
			sum += clusterWeights[i] * clusters[i].value(cPairs, beta);
		}
		return sum;
	}
	
	public ClusterAbstract[] cluster() {return clusters;}

	private final ClusterAbstract[] clusters;
    private final double[] clusterWeights;
}
