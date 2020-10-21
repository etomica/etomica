/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

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
		weightCoefficients = new double[allClusters.length];
		for (int i=0; i<weightCoefficients.length; i++) {
			weightCoefficients[i] = 1.0/weightCoefficients.length;
		}
	}
    
    public ClusterAbstract makeCopy() {
        ClusterAbstract[] copies = new ClusterAbstract[clusterArray.length];
        for (int i = 0; i < copies.length; i++) {
            copies[i] = clusterArray[i].makeCopy();
        }
        ClusterWeightUmbrella newCluster = new ClusterWeightUmbrella(copies);
        newCluster.setWeightCoefficients(weightCoefficients);
        return newCluster;
    }
	
	public int pointCount() {
		// can they be different?
		return clusterArray[0].pointCount();
	}

	public double value(BoxCluster box) {
		double sum = 0.0;
		for (int i=0; i<clusterArray.length; i++) {
			double v = clusterArray[i].value(box);
			sum += v*weightCoefficients[i];
		}
		return sum;
	}
	
	public void setWeightCoefficients(double[] a) {
	    System.arraycopy(a, 0, weightCoefficients, 0, a.length);
	}
	
	public double[] getWeightCoefficients() {
		return weightCoefficients;
	}
    
    public void setTemperature(double temp) {
        for (int i=0; i<clusterArray.length; i++) {
            clusterArray[i].setTemperature(temp);
        }
    }
    
    public ClusterAbstract[] getClusters() {
        return clusterArray;
    }
    
    private static final long serialVersionUID = 1L;
    private final ClusterAbstract[] clusterArray;
    private final double[] weightCoefficients;
}
