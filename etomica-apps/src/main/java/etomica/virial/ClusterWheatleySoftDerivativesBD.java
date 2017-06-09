/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial;

import etomica.math.SpecialFunctions;

import java.math.BigDecimal;
import java.math.MathContext;


/**
 * This class calculates the sum of all biconnected clusters using Wheatley's
 * recursive formulation.
 *
 * @author David Kofke and Andrew Schultz
 */
public class ClusterWheatleySoftDerivativesBD implements ClusterAbstract, ClusterAbstractMultivalue {

    protected final int n;
    protected final MayerFunction f;
    protected final int nDer;
    protected final BigDecimal[][] fQ, fC;
    protected final BigDecimal[][] fA, fB;
    public static BigDecimal BDONE = new BigDecimal(1.0);
    public static BigDecimal BDZERO = new BigDecimal(0.0);
    protected final MathContext mc;
    protected long cPairID = -1, lastCPairID = -1;
    protected final double[] value, lastValue;
    protected double beta;
    public static boolean pushme = false;
    protected double tol;
    protected ClusterWheatleySoftBD clusterBD;
    protected boolean debug = false;
    protected boolean doCaching = true;
    protected final BigDecimal[][] binomial;
    protected BigDecimal BDbeta;

    public ClusterWheatleySoftDerivativesBD(int nPoints, MayerFunction f, int precision, int nDer) {
        this.n = nPoints;
        value = new double[nDer + 1];
        lastValue = new double[nDer + 1];
        this.f = f;
        int nf = 1 << n;  // 2^n
        fQ = new BigDecimal[nf][nDer + 1];
        fC = new BigDecimal[nf][nDer + 1];
        for (int i = 0; i < n; i++) {
            for (int m = 0; m <= nDer; m++) {
                fQ[1 << i][m] = m == 0 ? BDONE : BDZERO;
            }
        }
        fA = new BigDecimal[nf][nDer + 1];
        fB = new BigDecimal[nf][nDer + 1];
        mc = new MathContext(precision);
        this.nDer = nDer;
        tol = Math.pow(10, -precision + 3);
        this.binomial = new BigDecimal[nDer + 1][];
        for (int m = 0; m <= nDer; m++) {
            binomial[m] = new BigDecimal[m + 1];
            for (int l = 0; l <= m; l++) {
                binomial[m][l] = new BigDecimal(SpecialFunctions.factorial(m) / (SpecialFunctions.factorial(l) * SpecialFunctions.factorial(m - l)));
            }
        }
    }

    public void setDoCaching(boolean newDoCaching) {
        doCaching = newDoCaching;
        if (clusterBD != null) {
            clusterBD.setDoCaching(doCaching);
        }
    }

    public ClusterAbstract makeCopy() {
        ClusterWheatleySoftDerivativesBD c = new ClusterWheatleySoftDerivativesBD(n, f, mc.getPrecision(), nDer);
        c.setTemperature(1 / beta);
        c.setDoCaching(doCaching);
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
                for (int m = 0; m <= nDer; m++) {
                    value[m] = lastValue[m];
                }
                return value[0];
            }

            // a new cluster
            lastCPairID = cPairID;
            for (int m = 0; m <= nDer; m++) {
                lastValue[m] = value[m];
            }
            cPairID = thisCPairID;
        }

        updateF(box);

        calcValue(box);
        if (Double.isNaN(value[0]) || Double.isInfinite(value[0])) {
            debug = true;
            updateF(box);
            calcValue(box);
            throw new RuntimeException("oops " + value[0]);
        }
        return value[0];
    }

    protected BigDecimal BDlog(BigDecimal val) {
        BigDecimal u = val.subtract(BDONE, mc);
        if (u.abs().doubleValue() > 0.01) {
            return new BigDecimal(Math.log(val.doubleValue()));
        } else {
            double x = -u.doubleValue();
            double tval = 0;
            double xn = 1;
            for (int i = 1; i < 11; i++) {
                xn *= x;
                tval -= xn / i;
            }
            return new BigDecimal(tval);
        }
    }

    /**
     * This calculates all FQ values given that the entries for pairs have
     * already been populated.
     */
    protected void calcFullFQ(BoxCluster box) {
        int nf = 1 << n;
        // generate all partitions and compute product of e-bonds for all pairs in partition
        for (int i = 3; i < nf; i++) {

            int j = i & -i;//lowest bit in i
            if (i == j) continue; // 1-point set
            int k = i & ~j; //strip j bit from i and set result to k
            if (k == (k & -k)) {
                if (fQ[i][0].doubleValue() == 0) {
                    for (int m = 1; m <= nDer; m++) {
                        fQ[i][m] = BDZERO;
                    }
                    continue;
                }
                BigDecimal c = BDlog(fQ[i][0]).divide(BDbeta, mc);
                for (int m = 1; m <= nDer; m++) {
                    fQ[i][m] = fQ[i][m - 1].multiply(c, mc);      //calculates derivatives of fQ w.r.t. beta
                }
                continue; // 2-point set; these fQ's were filled when bonds were computed, so skip
            }
            fQ[i][0] = fQ[k][0]; //initialize with previously-computed product of all pairs in partition, other than j
            if (fQ[i][0].doubleValue() == 0) {
                for (int m = 1; m <= nDer; m++) {
                    fQ[i][m] = BDZERO;
                }
                continue;
            }
            //loop over pairs formed from j and each point in partition; multiply by bond for each pair
            //all such pairs will be with bits higher than j, as j is the lowest bit in i
            for (int l = (j << 1); l < i; l = (l << 1)) {
                if ((l & i) == 0) continue; //l is not in partition
                fQ[i][0] = fQ[i][0].multiply(fQ[l | j][0], mc);
            }

            if (fQ[i][0].doubleValue() == 0) {
                for (int m = 1; m <= nDer; m++) {
                    fQ[i][m] = BDZERO;
                }
                continue;
            }
            BigDecimal c = BDlog(fQ[i][0]).divide(BDbeta, mc);
            for (int m = 1; m <= nDer; m++) {
                fQ[i][m] = fQ[i][m - 1].multiply(c, mc);      //calculates derivatives of fQ w.r.t. beta
            }
        }
    }

    /**
     * Returns the cluster value for the given configuration.  You must call
     * doCheck(BoxCluster) before calling this method.
     */
    public void calcValue(BoxCluster box) {
        double maxR2 = 0.1;
        if (pushme) {
            // force the system to hang out between minMaxR2 and maxMaxR2
            for (int i = 0; i < n - 1; i++) {
                for (int j = i + 1; j < n; j++) {
                    double r2 = box.getCPairSet().getr2(i, j);
                    if (r2 > maxR2) maxR2 = r2;
                }
            }
            double minMaxR2 = 5 * 5;
            if (maxR2 < minMaxR2) {
                value[0] = 1e-200;
                return;
            }
            double maxMaxR2 = 7 * 7;
            if (maxR2 > maxMaxR2) {
                value[0] = 0;
                return;
            }
        }
        calcFullFQ(box);
        int nf = 1 << n;
        //Compute the fC's
        for (int i = 1; i < nf; i++) {
            for (int m = 0; m <= nDer; m++) {
                fC[i][m] = fQ[i][m];
//        		System.out.println("fC["+i+"]["+m+"]"+fC[i][m]);
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

//                System.out.println("jcomp = "+jComp);
                for (int m = 0; m <= nDer; m++) {
                    for (int l = 0; l <= m; l++) {
                        fC[i][m] = fC[i][m].subtract(binomial[m][l].multiply(fC[j][l].multiply(fQ[jComp][m - l], mc), mc), mc);//for fQ, flip the bits on j; use only those appearing in i
                        //computes fC and its derivatives w.r.t beta
                    }
                }
            }
        }

        // find fA1
        for (int i = 2; i < nf; i += 2) {
            // all even sets don't contain 1
            //fA[i] = 0;
            for (int m = 0; m <= nDer; m++) {
                fB[i][m] = fC[i][m];
            }
        }
        for (int m = 0; m <= nDer; m++) {
            fA[1][m] = BDZERO;
            fB[1][m] = fC[1][m];
        }
        for (int i = 3; i < nf; i += 2) {
            // every set will contain 1
            for (int m = 0; m <= nDer; m++) {
                fA[i][m] = BDZERO;
                fB[i][m] = fC[i][m];
            }
            int ii = i - 1;//all bits in i but lowest
            int iLow2Bit = (ii & -ii);//next lowest bit
            int jBits = 1 | iLow2Bit;
            if (jBits == i) continue;
            //jBits has 1 and next lowest bit in i
            int iii = ii ^ iLow2Bit;//i with 2 lowest bits off
            int jInc = (iii & -iii);//3rd lowest bit, also increment for j
            for (int j = jBits; j < i; j += jInc) {//sum over partitions of i containing jBits
                int jComp = (i & ~j); //subset of i complementing j
                while ((j | jComp) != i && j < i) {
                    int jHighBits = j ^ jBits;
                    int jlow = jHighBits & -jHighBits;
                    j += jlow;
                    jComp = (i & ~j);
                }
                if (j == i) break;
                for (int m = 0; m <= nDer; m++) {
                    for (int l = 0; l <= m; l++) {
                        fA[i][m] = fA[i][m].add(binomial[m][l].multiply(fB[j][l].multiply(fC[jComp | 1][m - l], mc), mc), mc);
                    }
                }
            }
            for (int m = 0; m <= nDer; m++) {
                fB[i][m] = fB[i][m].subtract(fA[i][m], mc);//remove from B graphs that contain articulation point at 0
            }
        }

        for (int v = 1; v < n; v++) {
            int vs1 = 1 << v;
            for (int i = vs1 + 1; i < nf; i++) {
                for (int m = 0; m <= nDer; m++) {
                    fA[i][m] = BDZERO;
                }
//                fB[v][i] = fB[v-1][i];//no a.p. at v or below, starts with those having no a.p. at v-1 or below
                //rest of this is to generate A (diagrams having a.p. at v but not below), and subtract it from B
                if ((i & vs1) == 0) continue;//if i doesn't contain v, fA and fB are done
                int iLowBit = (i & -i);//lowest bit in i
                if (iLowBit == i) { //lowest bit is only bit; fA and fB are done
                    continue;
                }
                int jBits;
                int ii = i ^ iLowBit;
                int iLow2Bit = (ii & -ii);
                if (iLowBit != vs1 && iLow2Bit != vs1) {
                    //v is not in the lowest 2 bits
                    // jBits is the lowest bit and v
                    jBits = iLowBit | vs1;

                    // we can only increment by the 2nd lowest
                    int jInc = iLow2Bit;

                    //at this point jBits has (lowest bit + v) or (v + next lowest bit)
                    for (int j = jBits; j < i; j += jInc) {//sum over partitions of i
                        if ((j & jBits) != jBits) {
                            //ensure jBits are in j
                            j |= vs1;
                            if (j == i) break;
                        }
                        int jComp = i & ~j;//subset of i complementing j
                        while ((j | jComp) != i && j < i) {
                            int jHighBits = j ^ jBits;
                            int jlow = jHighBits & -jHighBits;
                            j += jlow; // this might knock out the v bit
                            j |= vs1;
                            jComp = (i & ~j);
                        }
                        if (j == i) break;
                        for (int m = 0; m <= nDer; m++) {
                            for (int l = 0; l <= m; l++) {
                                fA[i][m] = fA[i][m].add(binomial[m][l].multiply(fB[j][l].multiply((fB[jComp | vs1][m - l].add(fA[jComp | vs1][m - l], mc)), mc), mc), mc);
                            }
                        }
                    }
                } else {
                    //lowest 2 bits contain v
                    // jBits is the lowest 2 bits
                    // we can start at jBits and increment by the 3rd lowest bit
                    jBits = iLowBit | iLow2Bit;
                    if (jBits == i) continue; // no bits left for jComp

                    int iii = ii ^ iLow2Bit;
                    int jInc = (iii & -iii);

                    //at this point jBits has (lowest bit + v) or (v + next lowest bit)
                    for (int j = jBits; j < i; j += jInc) {//sum over partitions of i
                        // start=jBits and jInc ensure that every set includes jBits
                        int jComp = i & ~j;//subset of i complementing j
                        while ((j | jComp) != i && j < i) {
                            int jHighBits = j ^ jBits;
                            int jlow = jHighBits & -jHighBits;
                            j += jlow;
                            jComp = (i & ~j);
                        }
                        if (j == i) break;
                        for (int m = 0; m <= nDer; m++) {
                            for (int l = 0; l <= m; l++) {
                                fA[i][m] = fA[i][m].add(binomial[m][l].multiply(fB[j][l].multiply((fB[jComp | vs1][m - l].add(fA[jComp | vs1][m - l], mc)), mc), mc), mc);
                            }
                        }
                    }
                }
                for (int m = 0; m <= nDer; m++) {
                    fB[i][m] = fB[i][m].subtract(fA[i][m], mc);//remove from B graphs that contain articulation point at v
                }
            }
        }
        if (Math.abs(fB[nf - 1][0].doubleValue()) < tol) {
            if (clusterBD != null) {
                value[0] = clusterBD.value(box);
            } else {
                value[0] = 0;
            }
            return;
        }
//        System.out.println("fQ"+" "+Arrays.toString(fQ[nf-1]));
//        for(int i=0;i<fC.length;i++){
//            System.out.println("fC"+" "+i+" "+Arrays.toString(fC[i]));}
//        System.out.println("fA"+" "+Arrays.toString(fA[nf-1]));
//        System.out.println("fB"+" "+Arrays.toString(fB[nf-1]));
//        if (fB[nf-1][0] != -1 && fB[nf-1][1] == 0) {
//            throw new RuntimeException("oops");
//        }
        for (int m = 0; m <= nDer; m++) {
            value[m] = (1 - n) * fB[nf - 1][m].doubleValue() / SpecialFunctions.factorial(n);
        }
        if (pushme && maxR2 > 2 * 2) {
//            value *= Math.pow(maxR2/4, 6);
        }
    }

    protected void updateF(BoxCluster box) {
        CoordinatePairSet cPairs = box.getCPairSet();
        AtomPairSet aPairs = box.getAPairSet();

        f.setBox(box);
        // recalculate all f values for all pairs
        for (int i = 0; i < n - 1; i++) {
            for (int j = i + 1; j < n; j++) {
                double ff = f.f(aPairs.getAPair(i, j), cPairs.getr2(i, j), beta);
                if (false && Double.isNaN(ff)) {
                    f.f(aPairs.getAPair(i, j), cPairs.getr2(i, j), beta);
                    throw new RuntimeException("oops");
                }
//                if (Math.abs(ff) < 1e-14) ff = 0;
                fQ[(1 << i) | (1 << j)][0] = new BigDecimal(ff).add(BDONE, mc);
            }
        }

//        ANALYTICAL CHECK FOR B2 LJ
//        double sum = 0;
//        double sumpi = 0;
//        double dsum = 0;
//        double dr = 0.001;
//        for (int i=0; i<10000; i++) {
//            double r = dr*i;
//            double fval = f.f(null,  r*r, beta);
//            sum += fval*r*r*dr;
//            sumpi += Math.abs(fval)*r*r*dr;
//            double u = ((P2LennardJones)f.getPotential()).u(r*r);
//            if (u<Double.POSITIVE_INFINITY) {
//                dsum += -u*Math.exp(-u*beta)*r*r*dr;
//            }
//        }
//        System.out.println("B2 = "+(-0.5*4*Math.PI*sum) +" "+ (-0.5*4*Math.PI*dsum)+" "+(dsum/sumpi)+" "+2*Math.PI*sumpi);
    }

    public void setTemperature(double temperature) {
        beta = 1 / temperature;
        BDbeta = new BigDecimal(beta, mc);
        if (clusterBD != null) {
            clusterBD.setTemperature(temperature);
        }
    }

    public double[] getAllLastValues(BoxCluster box) {
        value(box);
        return value;
    }

    public static class ClusterRetrievePrimes implements ClusterAbstract {
        protected final ClusterWheatleySoftDerivativesBD cluster;
        protected final int n;

        public ClusterRetrievePrimes(ClusterWheatleySoftDerivativesBD cluster, int n) {
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
            return cluster.value[n];
        }

        public void setTemperature(double temperature) {


        }


    }

}
