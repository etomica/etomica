package etomica.virial;


public class ClusterSum implements ClusterAbstract {

	/**
	 * Constructor for ClusterSum.
	 */
	public ClusterSum(ClusterBonds[] subClusters, double[] subClusterWeights, MayerFunction[] fArray) {
        if (subClusterWeights.length != subClusters.length) throw new IllegalArgumentException("number of clusters and weights must be the same");
		clusters = new ClusterBonds[subClusters.length];
        clusterWeights = subClusterWeights;
        int pointCount = subClusters[0].pointCount();
		for(int i=0; i<clusters.length; i++) {
			clusters[i] = subClusters[i];
			if(clusters[i].pointCount() != pointCount) throw new IllegalArgumentException("Attempt to construct ClusterSum with clusters having differing numbers of points");
		}
        f = fArray;
        fValues = new double[pointCount][pointCount][fArray.length];
	}

	// equal point count enforced in constructor 
	public int pointCount() {
		return clusters[0].pointCount();
	}

	public double value(CoordinatePairSet cPairs, double beta) {
        int thisCPairID = cPairs.getID();
//        System.out.println(thisCPairID+" "+cPairID+" "+lastCPairID+" "+value+" "+lastValue+" "+f[0].getClass());
        if (thisCPairID == cPairID) return value;
        if (thisCPairID == lastCPairID) {
            // we went back to the previous cluster, presumably because the last
            // cluster was a trial that was rejected.  so drop the most recent value/ID
            cPairID = lastCPairID;
            value = lastValue;
            return value;
        }
        // a new cluster
        lastCPairID = cPairID;
        lastValue = value;
        cPairID = thisCPairID;
        
        int nPoints = pointCount();

        // recalculate all f values for all pairs
        for(int i=0; i<nPoints-1; i++) {
            for(int j=i+1; j<nPoints; j++) {
                for(int k=0; k<f.length; k++) {
                    fValues[i][j][k] = f[k].f(cPairs.getCPair(i,j),beta);
                    fValues[j][i][k] = fValues[i][j][k];
                }
            }
        }

		value = 0.0;
		for(int i=0; i<clusters.length; i++) {
            double v = clusters[i].value(fValues);
//            System.out.println("in cs.v "+clusterWeights[i]+" "+v);
            value += clusterWeights[i] * v;
//			value += clusterWeights[i] * clusters[i].value(fValues);
		}
//        System.out.println(value);
		return value;
	}
	
	public ClusterBonds[] cluster() {return clusters;}

	private final ClusterBonds[] clusters;
    private final double[] clusterWeights;
    private final MayerFunction[] f;
    private final double[][][] fValues;
    private int cPairID = -1, lastCPairID = -1;
    private double value, lastValue;
}