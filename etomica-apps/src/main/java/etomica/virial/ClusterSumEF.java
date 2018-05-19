/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial;

import java.util.Arrays;

/**
 * @author andrew
 *
 * TODO To change the template for this generated type comment go to
 * Window - Preferences - Java - Code Style - Code Templates
 */
public class ClusterSumEF extends ClusterSum {

    /**
     * Constructs a Cluster sum of the given cluster and weights.  The MayerFunction
     * array should contain only e-bonds.
     */
    public ClusterSumEF(ClusterBonds[] subClusters, double[] subClusterWeights, MayerFunction[] eArray) {
        super(subClusters, subClusterWeights, Arrays.copyOf(eArray, eArray.length * 2));
        numF = f.length/2;
    }

    public ClusterAbstract makeCopy() {
        MayerFunction[] e = new MayerFunction[numF];
        System.arraycopy(f,0,e,0,numF);
        ClusterSum copy = new ClusterSumEF(clusters,clusterWeights,e);
        copy.setTemperature(1/beta);
        return copy;
    }

    protected void updateF(BoxCluster box) {
        int nPoints = pointCount();
        for (int k=0; k<numF; k++) {
            f[k].setBox(box);
        }

        CoordinatePairSet cPairs = box.getCPairSet();
        AtomPairSet aPairs = box.getAPairSet();
        // recalculate all f values for all pairs
        for(int i=0; i<nPoints-1; i++) {
            for(int j=i+1; j<nPoints; j++) {
                for(int k=0; k<numF; k++) {
                    double eValue = f[k].f(aPairs.getAPair(i,j), cPairs.getr2(i,j), beta);
                    fValues[i][j][k+numF] = eValue;
                    fValues[j][i][k+numF] = eValue;
                    fValues[j][i][k] = eValue - 1;
                    fValues[i][j][k] = eValue - 1;
                }
            }
        }
    }

    private static final long serialVersionUID = 1L;
    private final int numF;
}
