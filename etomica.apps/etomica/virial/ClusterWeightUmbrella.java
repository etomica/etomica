package etomica.virial;


/**
 * Umbrella cluster weight wrapper function.  The umbrella cluster's
 * value will be positive when any of the sub clusters have non-zero values.
 * This differs from a ClusterSum mostly because the sub cluster values are
 * not actually summed.
 */
public class ClusterWeightUmbrella implements ClusterWeight {
	
	/**
	 * Contructs an umbrella cluster from the given clusters.
	 */
	public ClusterWeightUmbrella(ClusterSum[] allClusters) {
		clusterArray = allClusters;
		weightRatio = new double[allClusters.length];
		for (int i=0; i<weightRatio.length; i++) {
			weightRatio[i] = 1.0/weightRatio.length;
		}
	}
	
	public int pointCount() {
		// can they be different?
		return clusterArray[0].pointCount();
	}

	public double value(CoordinatePairSet cPairSet, AtomPairSet aPairSet, double beta) {
		double sum = 0.0;
		for (int i=0; i<clusterArray.length; i++) {
			double v = clusterArray[i].value(cPairSet,aPairSet,beta);
			sum += v*v*weightRatio[i];
		}
		return Math.sqrt(sum);
	}
	
	public void setWeightRatio(double[] aWeightRatio) {
		for (int i=weightRatio.length-1; i>=0; i--) {
			weightRatio[i] = aWeightRatio[i]/aWeightRatio[0];
		}
	}
	
	public double[] getWeightRatio() {
		return weightRatio;
	}
    
    private final ClusterAbstract[] clusterArray;
    private final double[] weightRatio;
}
