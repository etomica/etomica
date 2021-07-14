/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial.wheatley;

import etomica.math.SpecialFunctions;
import etomica.virial.AtomPairSet;
import etomica.virial.BoxCluster;
import etomica.virial.CoordinatePairSet;
import etomica.virial.MayerFunction;
import etomica.virial.cluster.ClusterAbstract;

import java.math.BigDecimal;
import java.math.MathContext;


/**
 * This class calculates the sum of all biconnected clusters using Wheatley's
 * recursive formulation.
 * 
 * @author David Kofke and Andrew Schultz 
 */
public class ClusterWheatleySoftBD implements ClusterAbstract {

    protected final int n;
    protected final MayerFunction f;

    protected final BigDecimal[] fQ, fC;
    protected final BigDecimal[] fA, fB;
    protected long cPairID = -1, lastCPairID = -1;
    protected double value, lastValue;
    protected double beta;
    public static boolean pushme = false;
    public static BigDecimal BDONE = new BigDecimal(1.0);
    public static BigDecimal BDZERO = new BigDecimal(0.0);
    protected final MathContext mc;
    protected double tol;
    protected boolean doCaching;

    public ClusterWheatleySoftBD(int nPoints, MayerFunction f, int precision) {
        this.n = nPoints;
        this.f = f;
        mc = new MathContext(precision);
        int nf = 1<<n;  // 2^n
        fQ = new BigDecimal[nf];
        fC = new BigDecimal[nf];
        for(int i=0; i<n; i++) {
            fQ[1<<i] = BDONE;
        }
        fA = new BigDecimal[nf];
        fB = new BigDecimal[nf];
        tol = Math.pow(10, -precision+3);
    }

    public void setDoCaching(boolean newDoCaching) {
        doCaching = newDoCaching;
    }

    public ClusterAbstract makeCopy() {
        ClusterWheatleySoftBD c = new ClusterWheatleySoftBD(n, f, mc.getPrecision());
        c.setTemperature(1/beta);
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
        }
      
        updateF(box);
      
        calcValue(box);
        if (Double.isNaN(value) || Double.isInfinite(value)) {
            updateF(box);
            calcValue(box);
            throw new RuntimeException("oops");
        }
        return value;
    }

    /**
     * This calculates all FQ values given that the entries for pairs have
     * already been populated.
     */
    protected void calcFullFQ(BoxCluster box) {
        int nf = 1<<n;
        // generate all partitions and compute product of e-bonds for all pairs in partition
        for (int i=3; i<nf; i++) {
            int j = i & -i;//lowest bit in i
            if (i==j) continue; // 1-point set
            int k = i&~j; //strip j bit from i and set result to k
            if (k == (k&-k)) continue; // 2-point set; these fQ's were filled when bonds were computed, so skip
            fQ[i] = fQ[k]; //initialize with previously-computed product of all pairs in partition, other than j
            if (fQ[i].doubleValue() == 0) continue;
            //loop over pairs formed from j and each point in partition; multiply by bond for each pair
            //all such pairs will be with bits higher than j, as j is the lowest bit in i
            for (int l=(j<<1); l<i; l=(l<<1)) {
                if ((l&i)==0) continue; //l is not in partition
                fQ[i] = fQ[i].multiply(fQ[l | j], mc);
            }
//            System.out.println(i+" "+fQ[i].toString());
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
            for (int i=0; i<n-1; i++) {
                for (int j=i+1; j<n; j++) {
                    double r2 = box.getCPairSet().getr2(i,j);
                    if (r2 > maxR2) maxR2 = r2;
                }
            }
            double minMaxR2 = 5*5;
            if (maxR2 < minMaxR2) {
                value = 1e-200;
                return;
            }
            double maxMaxR2 = 7*7;
            if (maxR2 > maxMaxR2) {
                value = 0;
                return;
            }
        }
        calcFullFQ(box);

        int nf = 1<<n;
        //Compute the fC's
//        System.out.println("fC");
        for(int i=1; i<nf; i++) {
            fC[i] = fQ[i];
            int iLowBit = i & -i;
            int inc = iLowBit<<1;
            for(int j=iLowBit; j<i; j+=inc) {
                int jComp = i & ~j;
                while ((j|jComp) != i && j<i) {
                    int jHighBits = j^iLowBit;
                    int jlow = jHighBits & -jHighBits;
                    j += jlow;
                    jComp = (i & ~j);
                }
                if (j==i) break;
                fC[i] = fC[i].subtract(fC[j].multiply(fQ[jComp], mc), mc);//for fQ, flip the bits on j; use only those appearing in i
            }
//            System.out.println(i+" "+fC[i].toString());
        }

        // find fA1
        for (int i=2; i<nf; i+=2) {
            // all even sets don't contain 1
            //fA[i] = 0;
            fB[i] = fC[i];
        }
        fA[1] = BDZERO;
        fB[1] = fC[1];
//        System.out.println("fA and fB");
        for (int i=3; i<nf; i+=2) {
            // every set will contain 1
            fA[i] = BDZERO;
            fB[i] = fC[i];
            int ii = i - 1;//all bits in i but lowest
            int iLow2Bit = (ii & -ii);//next lowest bit
            int jBits = 1 | iLow2Bit;
            if (jBits==i) continue;
            //jBits has 1 and next lowest bit in i
            int iii = ii ^ iLow2Bit;//i with 2 lowest bits off
            int jInc = (iii & -iii);//3rd lowest bit, also increment for j
            for (int j=jBits; j<i; j+=jInc) {//sum over partitions of i containing jBits
                int jComp = (i & ~j); //subset of i complementing j
                while ((j|jComp) != i && j<i) {
                    int jHighBits = j^jBits;
                    int jlow = jHighBits & -jHighBits;
                    j += jlow;
                    jComp = (i & ~j);
                }
                if (j==i) break;
                fA[i] = fA[i].add(fB[j].multiply(fC[jComp|1], mc), mc);
            }
            fB[i] = fB[i].subtract(fA[i], mc);//remove from B graphs that contain articulation point at 0
//            System.out.println("0 "+i+" "+fA[i]+" "+fB[i]);
        }

        for (int v=1; v<n; v++) {
            int vs1 = 1<<v;
            for (int i=vs1+1; i<nf; i++) {
                fA[i] = BDZERO;
//                fB[v][i] = fB[v-1][i];//no a.p. at v or below, starts with those having no a.p. at v-1 or below
                //rest of this is to generate A (diagrams having a.p. at v but not below), and subtract it from B
                if ((i & vs1) == 0) continue;//if i doesn't contain v, fA and fB are done
                int iLowBit = (i&-i);//lowest bit in i
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
                    for (int j=jBits; j<i; j+=jInc) {//sum over partitions of i
                        if ((j & jBits) != jBits) {
                            //ensure jBits are in j
                            j |= vs1;
                            if (j==i) break;
                        }
                        int jComp = i & ~j;//subset of i complementing j
                        while ((j|jComp) != i && j<i) {
                            int jHighBits = j^jBits;
                            int jlow = jHighBits & -jHighBits;
                            j += jlow; // this might knock out the v bit
                            j |= vs1;
                            jComp = (i & ~j);
                        }
                        if (j==i) break;
                        fA[i] = fA[i].add(fB[j].multiply(fB[jComp|vs1].add(fA[jComp|vs1], mc), mc), mc);
                    }
                }
                else {
                    //lowest 2 bits contain v
                    // jBits is the lowest 2 bits
                    // we can start at jBits and increment by the 3rd lowest bit
                    jBits = iLowBit | iLow2Bit;
                    if (jBits == i) continue; // no bits left for jComp
                    
                    int iii = ii ^ iLow2Bit;
                    int jInc = (iii & -iii);

                    //at this point jBits has (lowest bit + v) or (v + next lowest bit)
                    for (int j=jBits; j<i; j+=jInc) {//sum over partitions of i
                        // start=jBits and jInc ensure that every set includes jBits
                        int jComp = i & ~j;//subset of i complementing j
                        while ((j|jComp) != i && j<i) {
                            int jHighBits = j^jBits;
                            int jlow = jHighBits & -jHighBits;
                            j += jlow;
                            jComp = (i & ~j);
                        }
                        if (j==i) break;
                        fA[i] = fA[i].add(fB[j].multiply(fB[jComp|vs1].add(fA[jComp|vs1], mc), mc), mc);
                    }
                }

                fB[i] = fB[i].subtract(fA[i], mc);//remove from B graphs that contain articulation point at v
//                System.out.println(v+" "+i+" "+fA[i]+" "+fB[i]);
            }
        }
        if (Math.abs(fB[nf-1].doubleValue()) < tol) {
            value = 0;
            return;
        }
        value = (1-n)*fB[nf-1].doubleValue()/SpecialFunctions.factorial(n);
        if (pushme && maxR2 > 2*2) {
//            value *= Math.pow(maxR2/4, 6);
        }
    }

    protected void updateF(BoxCluster box) {
        CoordinatePairSet cPairs = box.getCPairSet();
        AtomPairSet aPairs = box.getAPairSet();

        f.setBox(box);
        // recalculate all f values for all pairs
        for(int i=0; i<n-1; i++) {
            for(int j=i+1; j<n; j++) {
                double ff = f.f(aPairs.getAPair(i,j),cPairs.getr2(i,j), beta);
//                if (Math.abs(ff) < 1e-14) ff = 0;
                fQ[(1<<i)|(1<<j)] = new BigDecimal(ff).add(BDONE, mc);
//                System.out.println(((1<<i)|(1<<j))+" "+ff+" "+fQ[(1<<i)|(1<<j)].toString());
            }
        }
    }

    public void setTemperature(double temperature) {
        beta = 1/temperature;
    }
}
