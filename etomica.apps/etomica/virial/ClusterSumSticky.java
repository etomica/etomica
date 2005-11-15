/*
 * Created on Oct 1, 2005
 *
 * TODO To change the template for this generated file go to
 * Window - Preferences - Java - Code Style - Code Templates
 */
package etomica.virial;

import etomica.simulation.Simulation;

/**
 * @author andrew
 *
 * TODO To change the template for this generated type comment go to
 * Window - Preferences - Java - Code Style - Code Templates
 */
public class ClusterSumSticky extends ClusterSum {

    public ClusterSumSticky(ClusterBonds[] subClusters, double[] subClusterWeights, MayerFunction[] fArray) {
        super(subClusters,subClusterWeights,fArray);
        iDiagram = 0; // full star, hopefully
        numDiagrams = clusters.length;
    }
    
    public double value(CoordinatePairSet cPairs, AtomPairSet aPairs) {
        if (isTrial) {
            //use cached f values, but recalculate the new diagram.
            isTrial = false;
            //forget the old lastValue since we will never go back to it.
            lastValue = value;
            lastCPairID = cPairID;
            calcValue();
            return value;
        }
        return super.value(cPairs,aPairs);
    }
    
    protected void calcValue() {
        value = numDiagrams * clusterWeights[iDiagram] * clusters[iDiagram].value(fValues);
    }        
    
    /**
     * set the diagram used for calculating the value to the ith diagram.
     */
    public void revertDiagram() {
//        System.out.println("switching back diagram "+iDiagram+" => "+oldDiagram);
        value = lastValue;
        cPairID = lastCPairID;
        iDiagram = oldDiagram;
    }

    /**
     * set the diagram used for calculating the value to a random one.
     */
    public void randomizeDiagram() {
        oldDiagram = iDiagram;
        iDiagram = Simulation.random.nextInt(numDiagrams);
        if (oldDiagram != iDiagram) {
            isTrial = true;
        }
//        System.out.println("switching from diagram "+oldDiagram+" => "+iDiagram);
    }

    private int iDiagram, oldDiagram;
    private boolean isTrial;
    private final int numDiagrams;
}
