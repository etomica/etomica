/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial;

import etomica.util.random.IRandom;

public class ClusterSumStickyEF extends ClusterSumEF {

    public ClusterSumStickyEF(ClusterBonds[] subClusters, double[] subClusterWeights,
            MayerFunction[] fArray, ClusterSumStickyEF[] buddies, IRandom random) {
        super(subClusters,subClusterWeights,fArray);
        this.random = random;
        iDiagram = 0; // full star, hopefully
        numDiagrams = clusters.length;
        clusterBuddies = buddies;
    }
    
    public ClusterSumStickyEF(ClusterSumEF cluster, ClusterSumStickyEF[] buddies, IRandom random) {
        this(cluster.clusters,cluster.clusterWeights,chopF(cluster.f),buddies, random);
    }
    
    public double value(BoxCluster box) {
        if (isTrial) {
            //use cached f values, but recalculate the new diagram.
            isTrial = false;
            //forget the old lastValue since we will never go back to it.
            lastValue = value;
            lastCPairID = cPairID;
            updateF(box);
            calcValue();
//            System.out.println("recalc diagram "+lastValue+" => "+value);
            return value;
        }
        return super.value(box);
    }
    
    protected void calcValue() {
        value = numDiagrams * clusterWeights[iDiagram] * clusters[iDiagram].value(fValues);
    }        
    
    /**
     * set the diagram used for calculating the value to the ith diagram.
     */
    public void revertDiagram() {
//        System.out.println("switching back diagram => "+oldDiagram); //+", value => "+lastValue);
        value = lastValue;
        cPairID = lastCPairID;
        iDiagram = oldDiagram;
        for (int i=0; i<clusterBuddies.length; i++) {
            clusterBuddies[i].revertDiagram();
        }
    }

    /**
     * set the diagram used for calculating the value to a random one.
     */
    public void randomizeDiagram() {
        oldDiagram = iDiagram;
        while (oldDiagram == iDiagram) {
            iDiagram = random.nextInt(numDiagrams);
        }
//        System.out.println("LJ => "+iDiagram);
        for (int i=0; i<clusterBuddies.length; i++) {
            clusterBuddies[i].setDiagram(iDiagram);
        }
        isTrial = true;
//        System.out.println("switching from diagram "+oldDiagram); //+" => "+iDiagram);
    }
    
    public void setDiagram(int i) {
//        System.out.println("hard sphere => "+i);
        oldDiagram = iDiagram;
        iDiagram = i;
        lastValue = value;
        lastCPairID = cPairID;
        isTrial = true;
    }

    static protected MayerFunction[] chopF(MayerFunction[] feArray) {
        MayerFunction[] eArray = new MayerFunction[feArray.length/2];
        System.arraycopy(feArray,0,eArray,0,eArray.length);
        return eArray;
    }
    
    private static final long serialVersionUID = 1L;
    public int iDiagram, oldDiagram;
    protected final IRandom random;
    private boolean isTrial;
    private final int numDiagrams;
    private final ClusterSumStickyEF[] clusterBuddies;
}
