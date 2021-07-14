/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial.cluster;

import etomica.virial.AtomPairSet;
import etomica.virial.BoxCluster;
import etomica.virial.CoordinatePairSet;
import etomica.virial.MayerFunction;

import java.util.Arrays;


/**
 * This class calculates the sum of all tree clusters using an adaptation of Wheatley's
 * recursive formulation.
 *
 * @author David Kofke and Andrew Schultz
 */
public class ClusterSinglyConnected implements ClusterAbstract {

    protected final int n, nf;
    protected final MayerFunction f;

    protected final double[] fL, fLN;
    protected final double[] bSum;
    protected long cPairID = -1, lastCPairID = -1;
    protected double value, lastValue;
    protected double beta;
    protected final int[][] partitions;

    public ClusterSinglyConnected(int nPoints, MayerFunction f) {
        this.n = nPoints;
        this.f = f;
        nf = 1<<n;  // 2^n
        fL = new double[nf];
        fLN = new double[nf];
        bSum = new double[nf];
        for (int i = 0; i < n; i++) {
            fL[1<<i] = 1;
            fLN[1<<i] = 1;
//            for(int j=i+1; j<n; j++) {
//                fNL[(1<<i)|(1<<j)] = 0;
//            }
        }
        partitions = new int[nf][];
        for (int m = 2; m < n; m++) {//structure as nested loops so we know what high bit is
            int iH = 1<<m; //high bit
            for (int i = iH + 3; i < (iH << 1); i++) {
                partitions[i] = computePartitions(i, iH);
            }
        }

    }

    protected final int[] computePartitions(int i, int iH) {

        int partitionCount = 0;
        int[] iPartitions = new int[10];//start with arbitrary length and resize as needed

        int iL = i & -i;//low bit
        int i0 = i^iL;//i, without the low bit
        if (i0 == iH) return null;//only two bits in i; we skip this because we start with all pairs in bSum and fL

        int inc = i0 & -i0;
        for (int iS = iH + iL; iS < i; iS += inc) {//structure loop to force iS to contain iH and iL bits
            int iSComp = i & ~iS;
            if ((iSComp | iS) != i) continue;
            partitionCount++;
            if (iPartitions.length < partitionCount) {
                iPartitions = Arrays.copyOf(iPartitions, partitionCount + 10);
            }
            iPartitions[partitionCount - 1] = iS;
        }

        if (iPartitions.length > partitionCount) {
            iPartitions = Arrays.copyOf(iPartitions, partitionCount);
        }
//        System.out.println(i+"\t"+Integer.bitCount(i)+"\t"+iPartitions.length);
        return iPartitions;
    }

    public ClusterAbstract makeCopy() {
        ClusterSinglyConnected c = new ClusterSinglyConnected(n, f);
        c.setTemperature(1 / beta);
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
                int index = (1 << i) | (1 << j);
                fLN[index] = fL[index] = bSum[index] = 1;
            }
        }

        calcValue();
        long num = Math.round(value);
        value = savedValue;
        return num;
    }

    /*
     * Computation of sum of tree diagrams.
     */
    protected void calcValue() {

        //Compute the fL and fN's
        // fL[i] is sum of all graphs in which low-bit node is a leaf
        // fNL[i] is sum of all graphs in which low-bit node is not a leaf or is a leaf
//        int sum1 = 0;
//        int sum[] = new int[n];
        for (int m = 2; m < n; m++) {//structure as nested loops so we know what high bit is
            final int iH = 1<<m; //high bit

            for (int i = iH + 3; i < (iH << 1); i++) {
                int iL = i & -i;//low bit
                int i0 = i^iL;//i, without the low bit
                if (i0 == iH) continue;//only two bits in i; we skip this because we start with all pairs in bSum and fL

                //bSum[i] is the sum of all bonds formed from the low-bit of i with each other non-zero bit of i.          
                //Calculation of bSum is performed for any i by adding the high-low (iH|iL) bit interaction to the sum
                //obtained without iH, computed from a previous iteration
                bSum[i] = bSum[i^iH] + bSum[iH|iL];
//                sum1++;
                //compute fNL and fL values
                fLN[i] = fL[i] = bSum[i] * fLN[i0];

                //add to fNL contributions from graphs where low bit is not a leaf
                int[] iPartitions = partitions[i];
                final int kmax = iPartitions.length;
                for (int k = 0; k < kmax; k++) {
                    int iS = iPartitions[k];
                    int iSComp = (i & ~iS);
                    fLN[i] += fL[iS] * fLN[iL | iSComp];
                    //sum[Integer.bitCount(i)-1]++;
                }
            }
        }

//        for(int k=0; k<n; k++) System.out.print("("+k+","+sum[k]+")\t");
//        System.out.println(sum1);
        value = fLN[nf-1];
    }

    protected void updateF(BoxCluster box) {
        CoordinatePairSet cPairs = box.getCPairSet();
        AtomPairSet aPairs = box.getAPairSet();

        f.setBox(box);
        // recalculate all f values for all pairs
        for (int i = 0; i < n - 1; i++) {
            for (int j = i + 1; j < n; j++) {
                int index = (1 << i) | (1 << j);
                bSum[index] = f.f(aPairs.getAPair(i, j), cPairs.getr2(i, j), beta);
                fLN[index] = fL[index] = bSum[index];
            }
        }
    }

    public void setTemperature(double temperature) {
        beta = 1 / temperature;
    }

    public static void main(String[] args) {
        for (int i = 0; i < 100; i++) {
            System.out.println(Integer.toBinaryString((1<<10)|i));
        }
        ClusterChainHS cc = new ClusterChainHS(5, null);
        ClusterSinglyConnected cs = new ClusterSinglyConnected(5, null);
        System.out.println(cc.numDiagrams());
        System.out.println(cs.numDiagrams());

    }
}
