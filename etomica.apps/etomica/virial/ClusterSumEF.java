/*
 * Created on Sep 24, 2005
 *
 * TODO To change the template for this generated file go to
 * Window - Preferences - Java - Code Style - Code Templates
 */
package etomica.virial;

import etomica.util.Arrays;

/**
 * @author andrew
 *
 * TODO To change the template for this generated type comment go to
 * Window - Preferences - Java - Code Style - Code Templates
 */
public class ClusterSumEF extends ClusterSum {

    /**
     * @param subClusters
     * @param subClusterWeights
     * @param fArray
     * @param temperature
     */
    public ClusterSumEF(ClusterBonds[] subClusters, double[] subClusterWeights,
            MayerFunction[] eArray, double temperature) {
        super(subClusters, subClusterWeights, (MayerFunction[])Arrays.resizeArray(eArray,eArray.length*2), temperature);
        numF = f.length/2;
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
    
    protected void updateF(CoordinatePairSet cPairs, AtomPairSet aPairs) {
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
                        eValue = ((MayerFunctionSpherical)f[k]).f(cPairs.getCPair(i,j),beta);
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
            return;
        }
        // recalculate all f values for all pairs
        for(int i=0; i<nPoints-1; i++) {
            for(int j=i+1; j<nPoints; j++) {
                for(int k=0; k<numF; k++) {
                    double eValue;
                    if (f[k] instanceof MayerFunctionSpherical) {
                        eValue = ((MayerFunctionSpherical)f[k]).f(cPairs.getCPair(i,j),beta);
                    }
                    else {
                        eValue = f[k].f(aPairs.getAPair(i,j),beta);
                    }
                    fValues[i][j][k+numF] = eValue;
                    fValues[j][i][k+numF] = eValue;
                    fValues[j][i][k] = eValue + 1;
                    fValues[i][j][k] = eValue + 1;
                }
            }
        }
    }
    
    private final int numF;
}
