/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial;

import etomica.math.SpecialFunctions;

/**
 * Cluster that uses recursion to compute temperature perturbation coefficients.
 * This awkwardly requires fR when it wants eR and fP when it wants uP.
 * This needs BD and some sort of criteria to determine when to use BD.
 * It acts a lot like ClusterWheatleyExtendSW, but seems more like derivatives
 * on paper.  It's missing the binomial coefficients that
 * ClusterWheatleySoftDerivatives uses.
 */
public class ClusterWheatleyPT implements ClusterAbstractMultivalue {

    protected long cPairID = -1, lastCPairID = -1;
    protected final double[] value, lastValue;
    protected final int n, nf, nOrder;
    protected final MayerFunction fR, fP;
    protected final double[][] fQ;
    protected final double[] fU;
    protected final double[][] fC, fA, fB;
    protected double beta;
    protected double tol;
    protected boolean doCaching = true;

    public ClusterWheatleyPT(int nPoints, int nOrder, MayerFunction fR, MayerFunction fP, double tol) {
        this.n = nPoints;
        this.nf = 1 << n; //2^n
        this.tol = tol;
        this.nOrder = nOrder;
        this.fR = fR;
        this.fP = fP;
        value = new double[nOrder + 1];
        lastValue = new double[nOrder + 1];
        fQ = new double[nf][nOrder + 1];
        for (int i = 0; i < n; i++) {
            fQ[1 << i][0] = 1.0;
            for (int if1 = 1; if1 <= nOrder; if1++) {
                fQ[1 << i][if1] = 0.0;
            }
        }
        fU = new double[nf];
        fC = new double[nf][nOrder + 1];
        fA = new double[nf][nOrder + 1];
        fB = new double[nf][nOrder + 1];
    }

    public ClusterAbstract makeCopy() {
        ClusterWheatleyPT c = new ClusterWheatleyPT(n, nOrder, fR, fP, tol);
        c.setTemperature(1 / beta);
        return c;
    }

    public int pointCount() {
        return n;
    }

    public double value(BoxCluster box) {
        if (doCaching) {
            CoordinatePairSet cPairs = box.getCPairSet();
            long thisCPairID = cPairs.getID();
            if (thisCPairID == cPairID) {
                return value[0];
            }
            if (thisCPairID == lastCPairID) {
                // we went back to the previous cluster, presumably because the last
                // cluster was a trial that was rejected.  so drop the most recent value/ID
                cPairID = lastCPairID;
                System.arraycopy(lastValue, 0, value, 0, value.length);
                return value[0];
            }

            // a new cluster
            lastCPairID = cPairID;
            System.arraycopy(value, 0, lastValue, 0, value.length);
            cPairID = thisCPairID;
        }
        updateF(box);

        calcValue(box);
        if (Double.isNaN(value[0]) || Double.isInfinite(value[0])) {
            updateF(box);
            calcValue(box);
            throw new RuntimeException("oops " + value[0]);
        }
        return value[0];
    }

    protected void updateF(BoxCluster box) {
        CoordinatePairSet cPairs = box.getCPairSet();
        AtomPairSet aPairs = box.getAPairSet();

        // recalculate all f values for all pairs
        for (int i = 0; i < n; i++) {
            for (int j = i + 1; j < n; j++) {
                int ij = 1 << i | 1 << j;
                fQ[ij][0] = fR.f(aPairs.getAPair(i, j), cPairs.getr2(i, j), beta) + 1;
                if (fQ[ij][0] == 0) {
                    fU[ij] = Double.POSITIVE_INFINITY;
                    for (int k = 1; k <= nOrder; k++) {
                        fQ[ij][k] = 0;
                    }
                    continue;
                }
                double fp = fP.f(aPairs.getAPair(i, j), cPairs.getr2(i, j), beta);

                double buij = -Math.log(fp + 1);
                fU[ij] = buij;
                for (int m = 1; m <= nOrder; m++) {
                    fQ[ij][m] = fQ[ij][m - 1] * (-buij) / m;
                }
            }
        }
    }

    /**
     * This calculates all FQ values given that the entries for pairs have
     * already been populated.
     */
    protected void calcFullFQ() {
        for (int i = 3; i < nf; i++) {
            int j = i & -i;
            if (i == j) continue; // 1-point set
            int k = i & ~j;
            if (k == (k & -k)) continue; // 2-point set

            fQ[i][0] = fQ[k][0];
            fU[i] = fU[k];

            //loop over pairs formed from j and each point in set i; multiply by bond for each pair
            //all such pairs will be with bits higher than j, as j is the lowest bit in i
            for (int l = (j << 1); l < i; l = (l << 1)) {
                if ((l & i) == 0) continue; //l is not in i
                fQ[i][0] = fQ[i][0] * fQ[l | j][0];
                fU[i] += fU[l | j];
            }
            if (fQ[i][0] == 0) {
                for (int m = 1; m <= nOrder; m++) {
                    fQ[i][m] = 0;
                }
                continue;
            }
            for (int m = 1; m <= nOrder; m++) {
                fQ[i][m] = fQ[i][m - 1] * (-fU[i] / m);
            }
        }
//		System.out.println("fQ="); MyUtility.display2DArray(fQ);
    }

    /**
     * Returns the cluster value for the given configuration.
     */
    public void calcValue(BoxCluster box) {
        calcFullFQ();
        //Compute the fC's
        for (int i = 1; i < nf; i++) {
            for (int if1 = 0; if1 <= nOrder; if1++) {
                fC[i][if1] = fQ[i][if1];
            }
            int iLowBit = i & -i;
            int inc = iLowBit << 1;

            for (int j = iLowBit; j < i; j += inc) {
                int jComp = i & ~j;
                while ((j | jComp) != i && j < i) {
                    int jHighBits = j ^ iLowBit;
                    int jlow = jHighBits & -jHighBits;
                    j += jlow;
                    jComp = (i & ~j);
                }
                if (j == i) break;
                for (int if1 = 0; if1 <= nOrder; if1++) {
                    for (int k = 0; k <= if1; k++) {
                        fC[i][if1] -= fC[j][k] * fQ[jComp][if1 - k];
                    }
                }
            }
        }
//       System.out.println("fC="); MyUtility.display2DArray(fC);

        //find fA1
        for (int i = 2; i < nf; i += 2) {//all even sets (2,4,6) don't contain 1
            for (int if1 = 0; if1 <= nOrder; if1++) {
                fB[i][if1] = fC[i][if1];
//    		   System.out.println("fB[" + i + "][" + if1+ "] =" + fB[i][if1]);
            }
        }

        for (int if1 = 0; if1 <= nOrder; if1++) {
            fA[1][if1] = 0;
            fB[1][if1] = fC[1][if1];
//           System.out.println("fB[" + 1 + "][" + if1+ "] =" + fB[1][if1]);
        }

        for (int i = 3; i < nf; i += 2) {//every set will contain 1.
            for (int if1 = 0; if1 <= nOrder; if1++) {
                fA[i][if1] = 0;
                fB[i][if1] = fC[i][if1];
            }
            int ii = i - 1;//all bits in i but lowest
            int iLow2Bit = (ii & -ii);//next lowest bit
            int jBits = 1 | iLow2Bit;
            if (jBits == i) continue;

            int iii = ii ^ iLow2Bit; //i with 2 lowest bits off
            int jInc = (iii & -iii);//3rd lowest bit, alsso increment for j
            for (int j = jBits; j < i; j += jInc) {//sum over partitions of i containing j Bits
                int jComp = (i & ~j);//subset of i complementing j
                while ((j | jComp) != i && j < i) {//if j is not a proper subset of i.
                    int jHighBits = j ^ jBits;
                    int jlow = jHighBits & -jHighBits;
                    j += jlow;
                    jComp = (i & ~j);
                }
                if (j == i) break;
//               if(i==nf-1) System.out.println("fA["+i+"]["+0+"]="+fA[i][0]+",fB["+i+"]["+0+"]="+fB[i][0]); 
                for (int if1 = 0; if1 <= nOrder; if1++) {
                    for (int k = 0; k <= if1; k++) {
                        fA[i][if1] += fB[j][k] * fC[jComp | 1][if1 - k];
//                     if((i==NF-1) && (if1==0)) printf("* %f, %f, %f\n", fA[i][if1], fB[j][k], fC[jComp|1][if1-k]);
//            		   if((i==nf-1) && (if1==0)) System.out.println("fA["+i+"]["+if1+"]="+fA[i][if1]+",[if1"+",fB["+j+"]["+k+"]="+fB[j][k]+",fC["+(jComp|1)+"]["+(if1-k)+"]="+fC[jComp|1][if1-k]);
                    }
                }
            }
            for (int if1 = 0; if1 <= nOrder; if1++) {
                fB[i][if1] -= fA[i][if1];//remove from B graphs that contain articulation point 0.
            }
//           if(i==nf-1) System.out.println("fB["+i+"]["+0+"]="+fB[i][0]); 
        }

        for (int v = 1; v < n; v++) {
            int vs1 = 1 << v;
            for (int i = vs1 + 1; i < nf; i++) {
                for (int if1 = 0; if1 <= nOrder; if1++) {
                    fA[i][if1] = 0;
                }
                if ((i & vs1) == 0) continue;
                int iLowBit = (i & -i);
                if (iLowBit == i) continue;

                int jBits;
                int ii = i ^ iLowBit;
                int iLow2Bit = (ii & -ii);
                if (iLowBit != vs1 && iLow2Bit != vs1) {
                    jBits = iLowBit | vs1;//v is not in the lowest 2 bits
                    int jInc = iLow2Bit;    //we can only increment by the 2nd lowest
                    for (int j = jBits; j < i; j += jInc) {
                        if ((j & jBits) != jBits) {
                            j |= vs1;
                            if (j == i) break;
                        }
                        int jComp = i & ~j;
                        while ((j | jComp) != i && j < i) {
                            int jHighBits = j ^ jBits;
                            int jlow = jHighBits & -jHighBits;
                            j += jlow;
                            j |= vs1;
                            jComp = (i & ~j);
                        }
                        if (j == i) break;
                        for (int if1 = 0; if1 <= nOrder; if1++) {
                            for (int k = 0; k <= if1; k++) {
                                fA[i][if1] += fB[j][k] * (fB[jComp | vs1][if1 - k] + fA[jComp | vs1][if1 - k]);
                            }
                        }
                    }
                } else {
                    //lowest 2 bits contain v
                    jBits = iLowBit | iLow2Bit;
                    if (jBits == i) continue; // no bits left jComp

                    int iii = ii ^ iLow2Bit;
                    int jInc = (iii & -iii);
                    //at this point jBits has (lowest bit + v)
                    for (int j = jBits; j < i; j += jInc) {//sum over partitions of i
                        int jComp = i & ~j;
                        while ((j | jComp) != i && j < i) {
                            int jHighBits = j ^ jBits;
                            int jlow = jHighBits & -jHighBits;
                            j += jlow;
                            jComp = (i & ~j);
                        }
                        if (j == i) break;
                        for (int if1 = 0; if1 <= nOrder; if1++) {
                            for (int k = 0; k <= if1; k++) {
                                fA[i][if1] += fB[j][k] * (fB[jComp | vs1][if1 - k] + fA[jComp | vs1][if1 - k]);
                            }
                        }
                    }
                }
                for (int if1 = 0; if1 <= nOrder; if1++) {
                    fB[i][if1] -= fA[i][if1];//remove from B graphs
                }
//        	   if(i==nf-1) System.out.println("fB["+i+"]["+0+"]=" + fB[i][0]);
            }
        }

        double fac = (1.0 - n) / SpecialFunctions.factorial(n);
        for (int i = 0; i < value.length; i++) {
            value[i] = fB[nf - 1][i] * fac;
        }
    }

    public void setTemperature(double temperature) {
        beta = 1 / temperature;
    }

    public ClusterRetrieveOrders makeClusterOrder(int o) {
        return new ClusterRetrieveOrders(o);
    }

    @Override
    public int getNumValues() {
        return nOrder + 1;
    }

    @Override
    public double[] getAllLastValues(BoxCluster box) {
        value(box);
        return value;
    }

    public class ClusterRetrieveOrders implements ClusterAbstract {
        protected final int iOrder;

        public ClusterRetrieveOrders(int iOrder) {
            this.iOrder = iOrder;
        }

        public ClusterAbstract makeCopy() {
            return null;
        }

        public int pointCount() {
            return n;
        }

        public double value(BoxCluster box) {
            return value[iOrder];
        }

        public void setTemperature(double temperature) {
        }
    }

    /**
     * Cluster returns only BD values.  configurations where BD was not used
     * result in value=0.
     */
    public static class ClusterRetrievePrimesBD implements ClusterAbstract {
        protected final ClusterWheatleySoftDerivatives cluster;
        protected final int n;

        public ClusterRetrievePrimesBD(ClusterWheatleySoftDerivatives cluster, int n) {
            this.cluster = cluster;
            this.n = n;
        }

        public ClusterAbstract makeCopy() {
            return null;
        }

        public int pointCount() {
            return cluster.n;
        }

        public double value(BoxCluster box) {
            return cluster.valueBD ? cluster.value[n] : 0;
        }

        public void setTemperature(double temperature) {
        }
    }
}
