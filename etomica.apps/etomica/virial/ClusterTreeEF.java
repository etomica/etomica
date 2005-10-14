/*
 * Created on Sep 24, 2005
 *
 * TODO To change the template for this generated file go to
 * Window - Preferences - Java - Code Style - Code Templates
 */
package etomica.virial;

import etomica.util.Arrays;
import etomica.virial.cluster.ClusterDiagramTree;

/**
 * @author andrew
 *
 * TODO To change the template for this generated type comment go to
 * Window - Preferences - Java - Code Style - Code Templates
 */
public class ClusterTreeEF extends ClusterTree {

    /**
     * @param bonds
     * @param fArray
     * @param temperature
     */
    public ClusterTreeEF(ClusterDiagramTree bonds, MayerFunction[] eArray,
            double temperature) {
        super(bonds, (MayerFunction[])Arrays.resizeArray(eArray,eArray.length*2), temperature);
        numF = f.length/2;
    }

    public ClusterAbstract makeCopy() {
        MayerFunction[] eArray = new MayerFunction[numF];
        System.arraycopy(f,0,eArray,0,numF);
        return new ClusterTreeEF(bondsTree,eArray,1/beta);
    }

    protected void updateF(CoordinatePairSet cPairs, AtomPairSet aPairs) {
        int nPoints = pointCount();
        
        if (cPairs.dirtyAtom > -1) {
            oldDirtyAtom = cPairs.dirtyAtom;
            int i = cPairs.dirtyAtom;
            int l = 0;
            for (int j=0; j<i; j++) {
                l += (i-j-1);
                for(int k=0; k<numF; k++) {
                    fOld[l][k] = fValues[l][k+numF];
                    if (f[k] instanceof MayerFunctionSpherical) {
                        fValues[l][k+numF] = ((MayerFunctionSpherical)f[k]).f(cPairs.getCPair(i,j),beta);
                    }
                    else {
                        fValues[l][k+numF] = f[k].f(aPairs.getAPair(i,j),beta);
                    }
                    fValues[l][k] = fValues[l][k] - 1;
                }
            }
            for (int j=i+1; j<nPoints; j++) {
                l++;
                for(int k=0; k<numF; k++) {
                    fOld[l][k] = fValues[l][k+numF];
                    if (f[k] instanceof MayerFunctionSpherical) {
                        fValues[l][k+numF] = ((MayerFunctionSpherical)f[k]).f(cPairs.getCPair(i,j),beta);
                    }
                    else {
                        fValues[l][k+numF] = f[k].f(aPairs.getAPair(i,j),beta);
                    }
                    fValues[l][k] = fValues[l][k+numF] - 1;
                }
            }
            return;
        }
        // recalculate all f values for all pairs
        int l=0;
        for(int i=0; i<nPoints-1; i++) {
            for(int j=i+1; j<nPoints; j++) {
                for(int k=0; k<numF; k++) {
                    if (f[k] instanceof MayerFunctionSpherical) {
                        fValues[l][k+numF] = ((MayerFunctionSpherical)f[k]).f(cPairs.getCPair(i,j),beta);
                    }
                    else {
                        fValues[l][k+numF] = f[k].f(aPairs.getAPair(i,j),beta);
                    }
                    fValues[l][k] = fValues[l][k+numF] - 1;
                }
                l++;
            }
        }
    }

    private final int numF;
}
