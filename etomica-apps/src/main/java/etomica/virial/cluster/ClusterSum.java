/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial.cluster;


import etomica.virial.AtomPairSet;
import etomica.virial.BoxCluster;
import etomica.virial.CoordinatePairSet;
import etomica.virial.MayerFunction;

import java.util.Arrays;

public class ClusterSum implements ClusterAbstract, java.io.Serializable {

    protected static final boolean debug = false;
    
    /**
     * Constructor for ClusterSum.  This class assumes that bonds defined in
     * subclusters with function indices that exceed the number of Mayer
     * functions are e-bonds.  With n Mayer functions, i=0..n-1 correspond
     * to the given Mayer functions (f), while i=n..2n-1 correspond to 
     * e=f+1, where f is the (i-n)th Mayer function.
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
        fOld = new double[pointCount][fArray.length];
        int maxF = 0;
        if (!clusters[0].isUsePermutations()) {
            // determine which fbonds are actually needed by the diagrams
            fullBondIndexArray = new int[pointCount-1][pointCount][0];
            for (int c=0; c<clusters.length; c++) {
                int[][] bondIndexArray = clusters[c].getBondIndexArray();
                for (int i=0; i<pointCount-1; i++) {
                    for (int j=i+1; j<pointCount; j++) {
                        int kf = bondIndexArray[i][j];
                        if (kf == -1) continue;
                        // if we have an ebond then first ensure that the fbond is listed; then return to ebond
                        int lmax = kf >= f.length ? 2 : 1;
                        kf = kf % f.length;
                        for (int l=0; l<lmax; l++) {
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
                                if (kf > maxF) maxF = kf;
                                fullBondIndexArray[i][j] = Arrays.copyOf(ff, ff.length + 1);
                                fullBondIndexArray[i][j][ff.length] = kf;
                            }
                            // we had an e-bond.  we need to remember to also calculate the f-bond
                            kf += f.length;
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
            for (int c=0; c<clusters.length; c++) {
                int[][] bondIndexArray = clusters[c].getBondIndexArray();
                for (int i=0; i<pointCount-1; i++) {
                    for (int j=i+1; j<pointCount; j++) {
                        int kf = bondIndexArray[i][j];
                        if (kf > maxF) maxF = kf;
                    }
                }
            }
        }
        fValues = new double[pointCount][pointCount][maxF+1];
    }
    
    public void setCaching(boolean doCaching) {
        this.doCaching = doCaching;
    }

    // equal point count enforced in constructor 
    public int pointCount() {
        return clusters[0].pointCount();
    }
    
    public ClusterAbstract makeCopy() {
        ClusterSum copy = new ClusterSum(clusters,clusterWeights,f);
        copy.setTemperature(1/beta);
        copy.setCaching(doCaching);
        return copy;
    }

    public double value(BoxCluster box) {
        if (doCaching) {
            CoordinatePairSet cPairs = box.getCPairSet();
            long thisCPairID = cPairs.getID();
         //  System.out.println(thisCPairID+" "+cPairID+" "+lastCPairID+" "+value+" "+lastValue+" "+f[0].getClass() +" here");
            if (thisCPairID == cPairID) {
            // System.out.println("clusterSum "+cPairID+" returning recent "+value);
                return value;
            }
            if (thisCPairID == lastCPairID) {
                // we went back to the previous cluster, presumably because the last
                // cluster was a trial that was rejected.  so drop the most recent value/ID
                cPairID = lastCPairID;
                value = lastValue;
             // System.out.println("clusterSum "+cPairID+" returning previous recent "+lastValue);
                return value;
            }

            // a new cluster
            lastCPairID = cPairID;
            lastValue = value;
            cPairID = thisCPairID;
        }
        
        updateF(box);
//        checkF(cPairs,aPairs);
        
        calcValue();
        if (debug && (Double.isNaN(value) || Double.isInfinite(value))) {
            updateF(box);
            calcValue();
            throw new RuntimeException("oops "+value);
        }
        return value;
    }
    
    protected void calcValue() {
        value = 0.0;
        for(int i=0; i<clusters.length; i++) {
            double v = clusters[i].value(fValues);
           // System.out.println(v + " v");
          //  if(v < 0){
                //System.out.println(" less than zero");
            //}
            value += clusterWeights[i] * v;
            System.out.println(v +" "+ clusterWeights[i] +" why v value " +  value +" value Product" );
          //  if(Double.isNaN(value)){
            //    System.out.println(debug + " requirement");
           // }
            // enable this to debug bogus values
            if (debug && (Double.isNaN(value) || Double.isInfinite(value)) || Double.isNaN(value)) {
                for (int j=0; j<fValues.length; j++) {
                    for (int k=0; k<fValues.length; k++) {
                      // System.out.println(j+" "+k+" "+Arrays.toString(fValues[j][k]));
                    }
                }
                throw new RuntimeException(value+" "+v);
            }
        }
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
                    if (fk < f.length) {
                        // we want the real fBond

                        fValues[i][j][fk] = f[fk].f(aPairs.getAPair(i,j),cPairs.getr2(i,j), beta);
                       // System.out.println(" ClusterSum Loop " +" " + fValues[i][j][fk]);
                    }
                    else {
                       // System.out.println("Entered else");
                        // we want an eBond
                        fValues[i][j][fk] = fValues[i][j][fk-f.length]+1;
                    }
                    if (debug && (Double.isNaN(fValues[i][j][fk]) || Double.isInfinite(fValues[i][j][fk]))) {
                       // System.out.println("entered f NAN");
                        f[fk].f(aPairs.getAPair(i,j),cPairs.getr2(i,j), beta);
                        throw new RuntimeException("oops f["+i+"]["+j+"]["+fk+"]="+fValues[i][j][fk]);
                    }
                    //System.out.println(fValues[i][j][fk] +" fvalues");
                    fValues[j][i][fk] = fValues[i][j][fk];
                }
            }
        }
       // System.out.println("out of the loop");
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
    protected int[][][] fullBondIndexArray;
    protected final MayerFunction[] f;
    protected double[][][] fValues;
    protected final double[][] fOld;
    protected long cPairID = -1, lastCPairID = -1;
    protected double value, lastValue;
    protected double beta;
    protected boolean doCaching = true;
}
