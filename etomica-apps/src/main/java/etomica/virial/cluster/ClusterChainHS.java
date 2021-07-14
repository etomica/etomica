/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial.cluster;


import etomica.virial.AtomPairSet;
import etomica.virial.BoxCluster;
import etomica.virial.CoordinatePairSet;
import etomica.virial.MayerFunction;

/**
 * This class calculates the sum of all chain and/or ring clusters for hard potentials, for which the Mayer function
 * can take values of only -1 or 0. Can return a general weighted sum of ring and chain values, special-casing
 * to chain-only and ring-only values if directed at construction.
 *
 * @author David Kofke and Andrew Schultz
 */
public class ClusterChainHS implements ClusterAbstract {

    protected final int n, nf1;
    protected final MayerFunction f;

    protected final double[][] nC;
    protected final double[][] fValues;
    protected double ringFrac, chainFrac;
    protected long cPairID = -1, lastCPairID = -1;
    protected double value, lastValue;
    public final boolean old = true;

    /**
     * Constructs with default to perform chain-only calculation, with chainFrac = 1.0
     */
    public ClusterChainHS(int nPoints, MayerFunction f) {
        this(nPoints, f, false);
    }

    /**
     * Constructs to perform either chain-only or ring-only calculation, as directed by doRing (true for ring-only)
     */
    public ClusterChainHS(int nPoints, MayerFunction f, boolean doRing) {
        this(nPoints, f, doRing?0:1, doRing?1:0);
    }

    /**
     * Constructs to perform linear combination of chain and ring values, with computed cluster values
     * each weighted by given fractions, and summed to get total value
     */
    public ClusterChainHS(int nPoints, MayerFunction f, double chainFrac, double ringFrac) {
        this.n = nPoints;
        this.f = f;
        this.chainFrac = chainFrac;
        this.ringFrac = ringFrac;
        nf1 = 1<<(n-1);  // 2^(n-1)

        nC = new double[n - 1][nf1];
        fValues = new double[n][n];
    }

    public ClusterAbstract makeCopy() {
        ClusterChainHS c = new ClusterChainHS(n, f);
        c.setTemperature(1);
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
        long num = Math.round(value);
        value = savedValue;
        return num;
    }

    /*
     * Computation of sum of chain or ring diagrams.
     */
    protected void calcValue() {

        //nC[m][i] is the number of chains beginning at last vertex and ending at m, traversing all points in i

        //Start with all pairwise paths from last vertex to each other vertex
        for(int m=0; m<n-1; m++) {
            nC[m][(1 << m)] = fValues[m][n - 1];
        }

        //All other paths
        for(int i=1; i<nf1; i++) {//i excludes the last bit, which is implicit in all partitions
            //the following two loops generate all pairs formed by each bit in i with each bit not in i

            //loop over bits not in i; start with full complement of i (i^(nf1-1)), and in each iteration
            //get lowest bit (im=(iC&-iC)) and strip it from complement (iC^=im) until complement is empty (iC=0)
            for(int iC=i^(nf1-1), im=(iC&-iC); iC>0; iC^=im,im=(iC&-iC)) {
                int m = log2(im);
                int iim = i|im;
                nC[m][iim] = 0;
             //loop over bits in i, in same manner as loop over complement
                for (int it = i, ik = (it & -it); ik > 0; it ^= ik, ik = (it & -it)) {
                    int k = log2(ik);
                    nC[m][iim] += fValues[m][k] * nC[k][i];
                }
            }//end for(iC)
        }//end for(i)

        double ringValue = 0;
        double chainValue = 0;

        if (ringFrac != 0.0) {
            for (int m = 0; m < n - 1; m++) {
                ringValue += nC[m][nf1-1] * fValues[m][n-1];
            }
        }

        if (chainFrac != 0.0) {

            //Sum chains in which last (n-1) vertex is not a leaf.
            //Consider all partitions, counting paths beginning in one partition and ending in its complement 
            //Use same looping structure as employed above
            for (int iS = 1; iS < nf1; iS += 2) {//keep 1 in iS-partition to prevent double counting
                //loop over bits not in iS
                int iSComp = iS^(nf1-1);
                for (int iC = iSComp, im = (iC & -iC); iC > 0; iC ^= im, im = (iC & -iC)) {
                    int m = log2(im);
                    //loop over bits in iS
                    for (int it = iS, ik = (it & -it); ik > 0; it ^= ik, ik = (it & -it)) {
                        int k = log2(ik);
                        chainValue += nC[m][iSComp] * nC[k][iS];
                    }
                }
            }

            //Sum chains where last (n-1) vertex is a leaf
            for (int m = 0; m < n - 1; m++) {
                chainValue += nC[m][nf1-1];
            }

        }//end if(chainFrac)

        value = chainFrac*chainValue + ringFrac*ringValue;
    }

    //gives position of bit for an integer having only one nonzero bit
    private final int log2(int i) {
        switch(i) {
            case (1 << 0):
                return 0;
            case (1 << 1):
                return 1;
            case (1 << 2):
                return 2;
            case (1 << 3):
                return 3;
            case (1 << 4):
                return 4;
            case (1 << 5):
                return 5;
            case (1 << 6):
                return 6;
            case (1 << 7):
                return 7;
            case (1 << 8):
                return 8;
            case (1 << 9):
                return 9;
            case (1 << 10):
                return 10;
            case (1 << 11):
                return 11;
            case (1 << 12):
                return 12;
            case (1 << 13):
                return 13;
            case (1 << 14):
                return 14;
            case (1 << 15):
                return 15;
            default:
                throw new IllegalArgumentException("Unexpected argument to log2: " + i);
        }
    }

    protected void updateF(BoxCluster box) {
        CoordinatePairSet cPairs = box.getCPairSet();
        AtomPairSet aPairs = box.getAPairSet();

        f.setBox(box);
        // recalculate all f values for all pairs
        for (int i = 0; i < n - 1; i++) {
            for (int j = i + 1; j < n; j++) {
                fValues[i][j] = f.f(aPairs.getAPair(i,j),cPairs.getr2(i,j), 1);
                fValues[j][i] = fValues[i][j];
            }
        }
    }

    public void setTemperature(double temperature) {
        // we don't need no stinkin temperature!
    }

    public static void main(String[] args) {

        for (int n = 5; n < 13; n++) {
//            ClusterChainWheatley cc = new ClusterChainWheatley(n, null);
            ClusterSinglyConnected cs = new ClusterSinglyConnected(n, null);
            ClusterChainHS cr = new ClusterChainHS(n, null, true);
            ClusterChainHS cc2 = new ClusterChainHS(n, null);
//            cc2.old = true;
            System.out.println(n + "\t" + cs.numDiagrams() +
                    "\t" + cr.numDiagrams() + "\t" + cc2.numDiagrams());
//            cc2.old = false;
            System.out.println(n + "\t" + cs.numDiagrams() +
                    "\t" + cr.numDiagrams() + "\t" + cc2.numDiagrams());

        }
    }
}
