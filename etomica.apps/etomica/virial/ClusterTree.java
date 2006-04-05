/*
 * Created on Sep 19, 2005
 *
 * TODO To change the template for this generated file go to
 * Window - Preferences - Java - Code Style - Code Templates
 */
package etomica.virial;

import etomica.virial.cluster.ClusterDiagramTree;

/**
 * @author andrew
 *
 * TODO To change the template for this generated type comment go to
 * Window - Preferences - Java - Code Style - Code Templates
 */
public class ClusterTree implements ClusterAbstract {

    public ClusterTree(ClusterDiagramTree bonds, MayerFunction[] fArray) {
        bondsTree = bonds;
        f = fArray;
        int nBody = bonds.pointCount();
        fValues = new double[nBody*(nBody-1)/2][fArray.length];
        fOld = new double[nBody][fArray.length];
    }

    public ClusterAbstract makeCopy() {
        ClusterTree copy = new ClusterTree(bondsTree,f);
        copy.setTemperature(1/beta);
        return copy;
    }

    public int pointCount() {
        return bondsTree.pointCount();
    }

    public double value(CoordinatePairSet cPairs, AtomPairSet aPairs) {
        int thisCPairID = cPairs.getID();
//      System.out.println(thisCPairID+" "+cPairID+" "+lastCPairID+" "+value+" "+lastValue+" "+f[0].getClass());
        if (thisCPairID == cPairID) {
            return value;
        }
        if (thisCPairID == lastCPairID) {
            // we went back to the previous cluster, presumably because the last
            // cluster was a trial that was rejected.  so drop the most recent value/ID
            if (oldDirtyAtom > -1) {
                System.out.println("reverting");
                revertF();
            }
            cPairID = lastCPairID;
            value = lastValue;
            return value;
        }
        // a new cluster
        lastCPairID = cPairID;
        lastValue = value;
        cPairID = thisCPairID;
        oldDirtyAtom = cPairs.dirtyAtom;
      
        updateF(cPairs,aPairs);
        value = bondsTree.value(fValues);
        return value;
    }
    
    protected void revertF() {
        int nPoints = pointCount();
        int i = oldDirtyAtom;
        int l = 0;
        for (int j=0; j<i; j++) {
            l += (i-j-1);
            for(int k=0; k<f.length; k++) {
                fValues[l][k] = fOld[j][k];
            }
            l += (nPoints-i);
        }
        for (int j=i+1; j<nPoints; j++) {
            for(int k=0; k<f.length; k++) {
                fValues[l][k] = fOld[j][k];
            }
            l++;
        }
        oldDirtyAtom = -1;
    }
    
    protected void updateF(CoordinatePairSet cPairs, AtomPairSet aPairs) {
        int nPoints = pointCount();
        
        if (cPairs.dirtyAtom > -1) {
            oldDirtyAtom = cPairs.dirtyAtom;
            int i = cPairs.dirtyAtom;
            int l = 0; //i-1;
            for (int j=0; j<i; j++) {
                l += (i-j-1);
                for(int k=0; k<f.length; k++) {
                    fOld[j][k] = fValues[l][k];
                    if (f[k] instanceof MayerFunctionSpherical) {
                        fValues[l][k] = ((MayerFunctionSpherical)f[k]).f(cPairs.getr2(i,j),beta);
                    }
                    else {
                        fValues[l][k] = f[k].f(aPairs.getAPair(i,j),beta);
                    }
                    if (Double.isInfinite(fValues[l][k])) {
                        System.out.println("oops9");
                    }
                }
                l += (nPoints-i);
            }
            for (int j=i+1; j<nPoints; j++) {
                for(int k=0; k<f.length; k++) {
                    fOld[j][k] = fValues[l][k];
                    if (f[k] instanceof MayerFunctionSpherical) {
                        fValues[l][k] = ((MayerFunctionSpherical)f[k]).f(cPairs.getr2(i,j),beta);
                    }
                    else {
                        fValues[l][k] = f[k].f(aPairs.getAPair(i,j),beta);
                    }
                    if (Double.isInfinite(fValues[l][k])) {
                        System.out.println("oops9");
                    }
                }
                l++;
            }
            return;
        }
        // recalculate all f values for all pairs
        int l=0;
        for(int i=0; i<nPoints-1; i++) {
            for(int j=i+1; j<nPoints; j++) {
                for(int k=0; k<f.length; k++) {
                    if (f[k] instanceof MayerFunctionSpherical) {
                        fValues[l][k] = ((MayerFunctionSpherical)f[k]).f(cPairs.getr2(i,j),beta);
                    }
                    else {
                        fValues[l][k] = f[k].f(aPairs.getAPair(i,j),beta);
                    }
                    if (Double.isInfinite(fValues[l][k])) {
                        System.out.println("oops9");
                    }
                }
                l++;
            }
        }
    }

    public void setTemperature(double temperature) {
        beta = 1/temperature;
    }
    
    public ClusterDiagramTree getBondsTree() {
        return bondsTree;
    }

    protected final ClusterDiagramTree bondsTree;
    protected final MayerFunction[] f;
    protected final double[][] fValues;
    private int cPairID = -1, lastCPairID = -1;
    private double value, lastValue;
    protected double beta;
    protected final double[][] fOld;
    protected int oldDirtyAtom;
}
