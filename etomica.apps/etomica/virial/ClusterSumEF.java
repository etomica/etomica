/*
 * Created on Sep 24, 2005
 *
 * TODO To change the template for this generated file go to
 * Window - Preferences - Java - Code Style - Code Templates
 */
package etomica.virial;

import etomica.atom.AtomPair;
import etomica.util.Arrays;

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
    public ClusterSumEF(ClusterBonds[] subClusters, double[] subClusterWeights,
            MayerFunction[] eArray) {
        super(subClusters, subClusterWeights, (MayerFunction[])Arrays.resizeArray(eArray,eArray.length*2));
        numF = f.length/2;
    }
    
    public ClusterAbstract makeCopy() {
        MayerFunction[] e = new MayerFunction[numF];
        System.arraycopy(f,0,e,0,numF);
        ClusterSum copy = new ClusterSumEF(clusters,clusterWeights,e);
        copy.setTemperature(1/beta);
        return copy;
    }

    protected void revertF() {
        int nPoints = pointCount();

        for(int j=0; j<nPoints; j++) {
            if (j == oldDirtyAtom) {
                continue;
            }
            for(int k=0; k<numF; k++) {
                double eValue = fOld[j][k];
                fValues[oldDirtyAtom][j][k+numF] = eValue;
                fValues[oldDirtyAtom][j][k] = eValue-1;
                fValues[j][oldDirtyAtom][k+numF] = eValue;
                fValues[j][oldDirtyAtom][k] = eValue-1;
            }
        }
    }
    
    protected void updateF(CoordinatePairLeafSet cPairs, AtomPairSet aPairs) {
        int nPoints = pointCount();

        if (cPairs.dirtyAtom > -1) {
            int i = cPairs.dirtyAtom;
            oldDirtyAtom = i;
            for(int j=0; j<nPoints; j++) {
                if (j == i) {
                    continue;
                }
                for(int k=0; k<numF; k++) {
                    double eValue;
                    //store the eValue
                    fOld[j][k] = fValues[i][j][k+numF];
                    if (f[k] instanceof MayerFunctionSpherical) {
                        double r2 = (i < j) ? cPairs.getr2(i,j) : cPairs.getr2(j,i);
                        eValue = ((MayerFunctionSpherical)f[k]).f(r2,beta);
                    }
                    else {
                        AtomPair pair = (i < j) ? aPairs.getAPair(i,j) : aPairs.getAPair(j,i);
                        eValue = f[k].f(pair,beta);
                    }
                    fValues[i][j][k+numF] = eValue;
                    fValues[j][i][k+numF] = eValue;
                    fValues[j][i][k] = eValue - 1;
                    fValues[i][j][k] = eValue - 1;
                }
            }
            return;
        }
        // recalculate all f values for all pairs
        for(int i=0; i<nPoints-1; i++) {
            for(int j=i+1; j<nPoints; j++) {
                for(int k=0; k<numF; k++) {
                    double eValue;
                    if (f[k] instanceof MayerFunctionSpherical) {
                        eValue = ((MayerFunctionSpherical)f[k]).f(cPairs.getr2(i,j),beta);
                    }
                    else {
                        eValue = f[k].f(aPairs.getAPair(i,j),beta);
                    }
                    fValues[i][j][k+numF] = eValue;
                    fValues[j][i][k+numF] = eValue;
                    fValues[j][i][k] = eValue - 1;
                    fValues[i][j][k] = eValue - 1;
                }
            }
        }
    }
    
    private final int numF;
}
