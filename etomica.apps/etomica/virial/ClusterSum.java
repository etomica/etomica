package etomica.virial;

import etomica.atom.AtomPair;



public class ClusterSum implements ClusterAbstract, java.io.Serializable {

    /**
     * Constructor for ClusterSum.
     */
    public ClusterSum(ClusterBonds[] subClusters, double[] subClusterWeights, MayerFunction[] fArray) {
        if (subClusterWeights.length != subClusters.length) throw new IllegalArgumentException("number of clusters and weights must be the same");
        clusters = new ClusterBonds[subClusters.length];
        clusterWeights = subClusterWeights;
        int pointCount = subClusters[0].pointCount();
        for(int i=0; i<clusters.length; i++) {
            clusters[i] = subClusters[i];
            if(clusters[i].pointCount() != pointCount) throw new IllegalArgumentException("Attempt to construct ClusterSum with clusters having differing numbers of points");
        }
        f = fArray;
        fValues = new double[pointCount][pointCount][fArray.length];
        fOld = new double[pointCount][fArray.length];
    }

    // equal point count enforced in constructor 
    public int pointCount() {
        return clusters[0].pointCount();
    }
    
    public ClusterAbstract makeCopy() {
        ClusterSum copy = new ClusterSum(clusters,clusterWeights,f);
        copy.setTemperature(1/beta);
        return copy;
    }

    public double value(CoordinatePairSet cPairs, AtomPairSet aPairs) {
        int thisCPairID = cPairs.getID();
//        System.out.println(thisCPairID+" "+cPairID+" "+lastCPairID+" "+value+" "+lastValue+" "+f[0].getClass());
        if (thisCPairID == cPairID) {
//            System.out.println("clusterSum returning recent "+value);
            return value;
        }
        if (thisCPairID == lastCPairID) {
            // we went back to the previous cluster, presumably because the last
            // cluster was a trial that was rejected.  so drop the most recent value/ID
            if (oldDirtyAtom > -1) {
                revertF();
            }
            cPairID = lastCPairID;
            value = lastValue;
//            System.out.println("clusterSum returning previous recent "+value);
            return value;
        }
        // a new cluster
        lastCPairID = cPairID;
        lastValue = value;
        cPairID = thisCPairID;
        
//        System.out.println("before");
        for (int i=0; i<pointCount()-1; i++) {
            for (int j=i+1; j<pointCount(); j++) {
//                System.out.println(i+" "+j+" "+fValues[i][j][0]);
            }
        }
        updateF(cPairs,aPairs);
//        System.out.println("updated");
        for (int i=0; i<pointCount()-1; i++) {
            for (int j=i+1; j<pointCount(); j++) {
//                System.out.println(i+" "+j+" "+fValues[i][j][0]);
            }
        }
//        checkF(cPairs,aPairs);
        
        calcValue();
        
//        System.out.println("clusterSum returning new "+value);
        return value;
    }
    
    protected void calcValue() {
        value = 0.0;
        for(int i=0; i<clusters.length; i++) {
            double v = clusters[i].value(fValues);
            value += clusterWeights[i] * v;
        }
    }        
    
    protected void revertF() {
        int nPoints = pointCount();

        for(int j=0; j<nPoints; j++) {
            if (j == oldDirtyAtom) {
                continue;
            }
            for(int k=0; k<f.length; k++) {
                fValues[oldDirtyAtom][j][k] = fOld[j][k];
                fValues[j][oldDirtyAtom][k] = fOld[j][k];
            }
        }
        oldDirtyAtom = -1;
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
                for(int k=0; k<f.length; k++) {
                    fOld[j][k] = fValues[i][j][k];
                    if (f[k] instanceof MayerFunctionSpherical) {
                        double r2 = (i < j) ? cPairs.getr2(i,j) : cPairs.getr2(j,i);
                        fValues[i][j][k] = ((MayerFunctionSpherical)f[k]).f(r2,beta);
                    }
                    else {
                        AtomPair pair = (i < j) ? aPairs.getAPair(i,j) : aPairs.getAPair(j,i);
                        fValues[i][j][k] = f[k].f(pair,beta);
                    }
                    fValues[j][i][k] = fValues[i][j][k];
                }
            }
            return;
        }
        // recalculate all f values for all pairs
        for(int i=0; i<nPoints-1; i++) {
            for(int j=i+1; j<nPoints; j++) {
                for(int k=0; k<f.length; k++) {
                    if (f[k] instanceof MayerFunctionSpherical) {
                        fValues[i][j][k] = ((MayerFunctionSpherical)f[k]).f(cPairs.getr2(i,j),beta);
                    }
                    else {
                        fValues[i][j][k] = f[k].f(aPairs.getAPair(i,j),beta);
                    }
                    fValues[j][i][k] = fValues[i][j][k];
                }
            }
        }
    }
    
    private void checkF(CoordinatePairSet cPairs, AtomPairSet aPairs) {
        int nPoints = pointCount();
        for(int i=0; i<nPoints-1; i++) {
            for(int j=i+1; j<nPoints; j++) {
                for(int k=0; k<f.length; k++) {
                    if (f[k] instanceof MayerFunctionSpherical) {
                        if (fValues[i][j][k] != ((MayerFunctionSpherical)f[k]).f(cPairs.getr2(i,j),beta)) {
                            throw new RuntimeException("oops1 "+i+" "+j+" "+k+" "+((MayerFunctionSpherical)f[k]).f(cPairs.getr2(i,j),beta));
                        }
                    }
                    else {
                        if (fValues[i][j][k] != f[k].f(aPairs.getAPair(i,j),beta)) {
                            throw new RuntimeException("oops2 "+i+" "+j+" "+k+" "+f[k].f(aPairs.getAPair(i,j),beta));
                        }
                    }
                    if (fValues[j][i][k] != fValues[i][j][k]) {
                        throw new RuntimeException("oops3 "+i+" "+j+" "+k+" "+fValues[j][i][k]+" "+fValues[i][j][k]);
                    }
                }
            }
        }
    }
    
    public ClusterBonds[] getClusters() {return clusters;}
    /**
     * @return Returns the temperature.
     */
    public double getTemperature() {
        return 1/beta;
    }
    /**
     * @param temperature The temperature to set.
     */
    public void setTemperature(double temperature) {
        beta = 1/temperature;
    }

    protected final ClusterBonds[] clusters;
    protected final double[] clusterWeights;
    protected final MayerFunction[] f;
    protected final double[][][] fValues;
    protected final double[][] fOld;
    protected int oldDirtyAtom;
    protected int cPairID = -1, lastCPairID = -1;
    protected double value, lastValue;
    protected double beta;
}
