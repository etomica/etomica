package etomica.virial;

import etomica.Atom;

/**
 * @author andrew
 *
 * cluster weight wrapper (absolute value of wrapped cluster)
 */
public class ClusterWeightAbs implements ClusterWeight {
	
	protected final ClusterAbstract weightCluster;
	
	public ClusterWeightAbs(ClusterAbstract cluster) {
		weightCluster = cluster;
	}
	
	public static ClusterWeight makeWeightCluster(ClusterAbstract cluster) {
		if (cluster instanceof ClusterWeight) {
			return (ClusterWeight)cluster;
		}
		return new ClusterWeightAbs(cluster);
	}

	public static ClusterWeight[] makeWeightClusters(ClusterAbstract[] clusters) {
		ClusterWeight[] weightClusters = new ClusterWeight[clusters.length];
		for (int i=0; i<clusters.length; i++) {
			weightClusters[i] = makeWeightCluster(clusters[i]);
		}
		return weightClusters;
	}
	
	public int pointCount() {
		return weightCluster.pointCount();
	}
	
	public double weight() {
		return weightCluster.weight();
	}
	
	public double value(CoordinatePairSet cPairSet, double beta) {
		return Math.abs(weightCluster.value(cPairSet, beta));
	}
	
}
