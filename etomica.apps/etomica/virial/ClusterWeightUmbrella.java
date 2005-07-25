package etomica.virial;


/**
 * Umbrella cluster weight wrapper function.  The umbrella cluster's
 * value will be positive when any of the sub clusters have non-zero values.
 * This differs from a ClusterSum mostly because the sub cluster values are
 * not actually summed.
 */
public class ClusterWeightUmbrella implements ClusterWeight, java.io.Serializable {
	
	/**
	 * Contructs an umbrella cluster from the given clusters.
	 */
	public ClusterWeightUmbrella(ClusterAbstract[] allClusters) {
		clusterArray = allClusters;
		weightRatio = new double[allClusters.length];
		for (int i=0; i<weightRatio.length; i++) {
			weightRatio[i] = 1.0/weightRatio.length;
		}
	}
    
    public ClusterAbstract makeCopy() {
        ClusterWeightUmbrella newCluster = new ClusterWeightUmbrella(clusterArray);
        newCluster.setWeightRatio(weightRatio);
        return newCluster;
    }
	
	public int pointCount() {
		// can they be different?
		return clusterArray[0].pointCount();
	}

	public double value(CoordinatePairSet cPairSet, AtomPairSet aPairSet) {
		double sum = 0.0;
		for (int i=0; i<clusterArray.length; i++) {
			double v = clusterArray[i].value(cPairSet,aPairSet);
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
    
    public void setTemperature(double temp) {
        for (int i=0; i<clusterArray.length; i++) {
            clusterArray[i].setTemperature(temp);
        }
    }        
    
    private final ClusterAbstract[] clusterArray;
    private final double[] weightRatio;
}
