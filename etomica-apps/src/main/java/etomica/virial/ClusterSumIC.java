/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial;

import etomica.box.Box;
import etomica.atom.IMoleculeList;
import etomica.api.IPotential;


/**
 * @author andrew
 *
 * cluster weight wrapper (absolute value of wrapped cluster)
 */
public class ClusterSumIC extends ClusterSum {
	
    private static final long serialVersionUID = 1L;
    protected final MayerFunction[] f12;
    protected final double[] subClusterValues;
    protected final double[] lastSubClusterValues;
	
    /**
     * @param subClusters array of clusters corresponding to distribution functions
     * @param subClusterWeights weights for the clusters
     * @param fArray Mayer functions used in the clusters.  the last one should be an
     *        instance of FR2.  each subCluster should have an FR2 bond between the
     *        first two points.
     * @param f12 array of Mayer functions (between the first two points) for each cluster 
     */
	public ClusterSumIC(ClusterBonds[] subClusters, double[] subClusterWeights, MayerFunction[] fArray, MayerFunction[] f12) {
	    super(subClusters, subClusterWeights, fArray);
		this.f12 = f12;
		subClusterValues = new double[subClusters.length+1];
        lastSubClusterValues = new double[subClusters.length+1];
	}
		
    public ClusterAbstract makeCopy() {
        ClusterSumIC copy = new ClusterSumIC(clusters, clusterWeights, f, f12);
        copy.setTemperature(1/beta);
        copy.setCaching(doCaching);
        return copy;
    }

    public int getNumValues() {
        return subClusterValues.length;
    }
    
    public double value12(int i) {
        return subClusterValues[i];
    }
    
    public double value(BoxCluster box) {
        CoordinatePairSet cPairs = box.getCPairSet();
        if (doCaching) {
            long thisCPairID = cPairs.getID();
//            System.out.println(thisCPairID+" "+cPairID+" "+lastCPairID+" "+value+" "+lastValue+" "+f[0].getClass());
            if (thisCPairID == cPairID) {
//                System.out.println("clusterSum "+cPairID+" returning recent "+value);
                return super.value(box);
            }
            if (thisCPairID == lastCPairID) {
                // we went back to the previous cluster, presumably because the last
                // cluster was a trial that was rejected.  so drop the most recent value/ID
                System.arraycopy(lastSubClusterValues, 0, subClusterValues, 0, subClusterValues.length);
//                System.out.println("clusterSum "+cPairID+" returning previous recent "+lastValue);
                return super.value(box);
            }

            // a new cluster
            System.arraycopy(subClusterValues, 0, lastSubClusterValues, 0, subClusterValues.length);
        }
        return super.value(box);
    }
    
    protected void calcValue() {
        value = 0.0;
        double r12 = fValues[0][1][f.length-1];
        if (r12==0) {
            value = 0;
            return;
        }
        subClusterValues[clusters.length] = 0;
        for(int i=0; i<clusters.length; i++) {
            double v = clusters[i].value(fValues)/r12;
            subClusterValues[i] = v;
            double x = clusterWeights[i] * v * f12[i].f(null, r12, beta);
            subClusterValues[clusters.length] += Math.abs(x);
            value += x;
        }
    }
    
    public static class FR2 implements MayerFunction {

        public double f(IMoleculeList pair, double r2, double beta) {
            return r2;
        }

        public IPotential getPotential() {return null;}
        public void setBox(Box box) {}
    }
}
