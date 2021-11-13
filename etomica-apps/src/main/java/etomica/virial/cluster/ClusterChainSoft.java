/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial.cluster;


import etomica.virial.AtomPairSet;
import etomica.virial.BoxCluster;
import etomica.virial.CoordinatePairSet;
import etomica.virial.MayerFunction;

/**
 * This class calculates the sum of all chain clusters. Can be configured at construction to do rings instead.
 * 
 * @author David Kofke and Andrew Schultz 
 */
public class ClusterChainSoft implements ClusterAbstract {

    protected final int n, nf;
    protected final MayerFunction f;
    protected final boolean doRing;
    
    protected final double[][] nC;
    protected final double[][] fValues;
    protected long cPairID = -1, lastCPairID = -1;
    protected double value, lastValue;
    protected double beta;
    
    public ClusterChainSoft(int nPoints, MayerFunction f) {
        this(nPoints, f, false);
    }
    
    public ClusterChainSoft(int nPoints, MayerFunction f, boolean doRing) {
        this.n = nPoints;
        this.f = f;
        this.doRing = doRing;
        nf = 1<<n;  // 2^n

        nC = new double[n][nf];
        fValues = new double[n][n];
    }

    public ClusterAbstract makeCopy() {
        ClusterChainSoft c = new ClusterChainSoft(n, f);
        c.setTemperature(1/beta);
        return c;
    }
    
    public int pointCount() {
        return n;
    }


    public double value(BoxCluster box) {
        CoordinatePairSet cPairs = box.getCPairSet();
        long thisCPairID = cPairs.getID();
        if (thisCPairID == cPairID) {
            return value;
        }
        if (thisCPairID == lastCPairID) {
            // we went back to the previous cluster, presumably because the last
            // cluster was a trial that was rejected.  so drop the most recent value/ID
            cPairID = lastCPairID;
            value = lastValue;
            return value;
        }

        // a new cluster
        lastCPairID = cPairID;
        lastValue = value;
        cPairID = thisCPairID;
        
        updateF(box);
        
        calcValue();
        if (Double.isNaN(value) || Double.isInfinite(value)) {
            updateF(box);
            calcValue();
            throw new RuntimeException("oops");
        }
        return value;
    }

    public long numDiagrams() {
        double savedValue = value;
        for (int i=0; i<n; i++) {
            for (int j=i+1; j<n; j++) {
                 fValues[i][j] = 1;
                 fValues[j][i] = 1;
            }
        }
        
        calcValue();
        long num = (int)Math.round(value);
        value = savedValue;
        return num;
    }


    /*
     * Computation of sum of chain or ring diagrams.
     */
    protected void calcValue() {
        
        //nC[m][i] is the number of chains beginning at 0 and ending at m, traversing all points in i
        
        //Start with all pairwise paths from 0 to each vertex
        for(int m=1; m<n; m++) {
            nC[m][(1<<m)|1] = fValues[m][0];
        }
        
        //All other paths
        //(could probably reduce memory by not including redundant first bit)
        for(int i=3; i<nf-1; i+=2) {//1-bit is always nonzero in i
            for(int m=1; m<n; m++) {//loop over indices not in i
                int im = 1<<m;
                if((im & i) != 0) continue;//skip if m is in i
                int index = i|im;
                nC[m][index] = 0;
                for(int k=1; k<n; k++) {//loop over indices in i
                    int ik = 1<<k;
                    if(ik > i) break;
                    if((ik & i) == 0) continue;//skip if k is not in i
                    nC[m][index] += fValues[m][k]*nC[k][i];
                }
            }
        }
        
        value = 0;
        if(doRing) {
            for(int m=1; m<n; m++) {
                value += nC[m][nf-1] * fValues[m][0];
            }
            
        } else { //chains
        
            //Sum chains in which first vertex is not a leaf.
            //Consider all partitions, counting paths beginning in one partition and ending in its complement
            for(int iS=3; iS<nf; iS+=4) {//keep 1 and 2 in i-partition to prevent double counting
                int iSComp = (nf-1)^iS;
                for(int m=1; m<n; m++) {
                    if(((1<<m)&iS) == 0) continue;//skip if m is not in iS
                    for(int k=2; k<n; k++) {
                        if(((1<<k)&iSComp) == 0) continue;//skip if k is not in iSComp
                        value += nC[m][iS] * nC[k][iSComp|1];
                    }
                }
            }
            
            //Sum chains where first vertex is a leaf
            for(int m=1; m<n; m++) {
                value += nC[m][nf-1];
            }
        
        }
        
    }
        
    protected void updateF(BoxCluster box) {
        CoordinatePairSet cPairs = box.getCPairSet();
        AtomPairSet aPairs = box.getAPairSet();

        f.setBox(box);
        // recalculate all f values for all pairs
        for(int i=0; i<n-1; i++) {
            for(int j=i+1; j<n; j++) {
                fValues[i][j] = f.f(aPairs.getAPair(i,j),cPairs.getr2(i,j), beta);
                fValues[j][i] = fValues[i][j];
            }
        }
    }
        
    public void setTemperature(double temperature) {
        beta = 1/temperature;
    }
}
