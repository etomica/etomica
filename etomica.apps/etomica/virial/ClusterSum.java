package etomica.virial;

import etomica.util.Arrays;
import etomica.util.Debug;


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
        
        if (!clusters[0].isUsePermutations()) {
            // determine which fbonds are actually needed by the diagrams
            fullBondIndexArray = new int[pointCount-1][pointCount][0];
            for (int c=0; c<clusters.length; c++) {
                int[][] bondIndexArray = clusters[c].getBondIndexArray();
                for (int i=0; i<pointCount-1; i++) {
                    for (int j=i+1; j<pointCount; j++) {
                        int kf = bondIndexArray[i][j];
                        if (kf == -1) continue;
                        int[] ff = fullBondIndexArray[i][j];
                        boolean newF = true;
                        for (int k=0; k<ff.length; k++) {
                            if (ff[k] == kf) {
                                // we'll already calculate MayerFunction kf for the i-j pair
                                newF = false;
                                break;
                            }
                        }
                        if (newF) {
                            fullBondIndexArray[i][j] = Arrays.resizeArray(ff, ff.length+1);
                            fullBondIndexArray[i][j][ff.length] = kf;
                        }
                    }
                }
            }
        }
        else {
            // when using permutations in ClusterBonds, everything will get rearranged
            // at some point, so each pair will have each bond
            fullBondIndexArray = new int[pointCount-1][pointCount][f.length];
            for (int i=0; i<pointCount-1; i++) {
                for (int j=i+1; j<pointCount; j++) {
                    for (int k=0; k<f.length; k++) {
                        fullBondIndexArray[i][j][k] = k;
                    }
                }
            }
        }
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

    public double value(BoxCluster box) {
        CoordinatePairSet cPairs = box.getCPairSet();
        int thisCPairID = cPairs.getID();
//        System.out.println(thisCPairID+" "+cPairID+" "+lastCPairID+" "+value+" "+lastValue+" "+f[0].getClass());
        if (thisCPairID == cPairID) {
//            System.out.println("clusterSum "+cPairID+" returning recent "+value);
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
//            System.out.println("clusterSum "+cPairID+" returning previous recent "+lastValue);
            return value;
        }

        // a new cluster
        lastCPairID = cPairID;
        lastValue = value;
        cPairID = thisCPairID;
        
        updateF(box);
//        checkF(cPairs,aPairs);
        
        calcValue();
//        if (Double.isNaN(value) || Double.isInfinite(value)) {
//            throw new RuntimeException("oops");
//        }
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
    
    protected void updateF(BoxCluster box) {
        int nPoints = pointCount();
        CoordinatePairSet cPairs = box.getCPairSet();
        AtomPairSet aPairs = box.getAPairSet();

        for (int k=0; k<f.length; k++) {
            f[k].setBox(box);
        }
        // recalculate all f values for all pairs
        for(int i=0; i<nPoints-1; i++) {
            for(int j=i+1; j<nPoints; j++) {
                // only update the mayer functions that we'll need for this pair
                int[] fij = fullBondIndexArray[i][j];
                for(int k=0; k<fij.length; k++) {
                    int fk = fij[k];
                    fValues[i][j][fk] = f[fk].f(aPairs.getAPair(i,j),cPairs.getr2(i,j), beta);
                    fValues[j][i][fk] = fValues[i][j][fk];
                }
            }
        }
    }
    
    private void checkF(BoxCluster box) {
        CoordinatePairSet cPairs = box.getCPairSet();
        AtomPairSet aPairs = box.getAPairSet();
        int nPoints = pointCount();
        for(int i=0; i<nPoints-1; i++) {
            for(int j=i+1; j<nPoints; j++) {
                int[] fij = fullBondIndexArray[i][j];
                for(int k=0; k<fij.length; k++) {
                    int fk = fij[k];
                    if (fValues[i][j][fk] != f[fk].f(aPairs.getAPair(i,j),cPairs.getr2(i,j), beta)) {
                        throw new RuntimeException("oops2 "+i+" "+j+" "+fk+" "+f[fk].f(aPairs.getAPair(i,j),cPairs.getr2(i,j), beta));
                    }
                    if (fValues[j][i][fk] != fValues[i][j][fk]) {
                        throw new RuntimeException("oops3 "+i+" "+j+" "+fk+" "+fValues[j][i][fk]+" "+fValues[i][j][fk]);
                    }
                }
            }
        }
    }
    
    public ClusterBonds[] getClusters() {return clusters;}
    public double[] getWeights() {return clusterWeights;}
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

    public double[][][] getFValues() {
        return fValues;
    }

    private static final long serialVersionUID = 1L;
    protected final ClusterBonds[] clusters;
    protected final double[] clusterWeights;
    int[][][] fullBondIndexArray;
    protected final MayerFunction[] f;
    protected double[][][] fValues;
    protected final double[][] fOld;
    protected int oldDirtyAtom;
    protected int cPairID = -1, lastCPairID = -1;
    protected double value, lastValue;
    protected double beta;
}
