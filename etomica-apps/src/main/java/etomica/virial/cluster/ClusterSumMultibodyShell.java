/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial.cluster;

import etomica.virial.BoxCluster;
import etomica.virial.MayerFunction;
import etomica.virial.MayerFunctionNonAdditive;

public class ClusterSumMultibodyShell extends ClusterSumMultibody {

    /**
     * Constructs a cluster (shell) that depends on another ClusterSumMultibody instance
     * (the core) to do the actual work of calculating all the f values.  This
     * class still does the work of multiplying f bonds together and summing
     * diagrams together.
     * 
     * MayerFunction arrays are needed only to make the ClusterSumMultibody superclass
     * happy.  Length of array should probably match core's f array.
     */
    public ClusterSumMultibodyShell(ClusterSumMultibody coreCluster, ClusterBonds[] subClusters,
                                    double[] subClusterWeights, MayerFunction[] f, MayerFunctionNonAdditive[] fMulti) {
        super(subClusters, subClusterWeights, f, fMulti);
        // grab this now and we won't need to do this again
        fValues = coreCluster.getFValues();
        fNonAdditiveValues = coreCluster.getFNonAdditiveValues();
        this.coreCluster = coreCluster;
    }
    
    public ClusterAbstract makeCopy() {
        throw new RuntimeException("this will only lead to pain");
    }
    
    protected void updateF(BoxCluster box) {
        // our f values are always up to date
    }
    
    protected final ClusterSumMultibody coreCluster;
}
