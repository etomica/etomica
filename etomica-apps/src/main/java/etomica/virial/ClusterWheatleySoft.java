/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial;

import etomica.math.SpecialFunctions;
import etomica.util.random.IRandom;


/**
 * This class calculates the sum of all biconnected clusters using Wheatley's
 * recursive formulation.
 * 
 * @author David Kofke and Andrew Schultz 
 */
public class ClusterWheatleySoft implements ClusterAbstract {

    protected final int n;
    protected final MayerFunction f;

    protected final double[] fQ, fC;
    protected final double[] fA, fB;
    protected long cPairID = -1, lastCPairID = -1;
    protected double value, lastValue;
    protected double beta;
    public static boolean pushme = false;
    protected double tol;
    protected ClusterWheatleySoftBD clusterBD;
    protected boolean debug = false;
    protected boolean doCaching = true;
    protected long totCount, countBD;
    protected double BDAccFrac;
    protected long timeBD;
    protected IRandom random;
    protected double[] avgAbsCheck = {0, 0}, avgAbsCheckBD = {0, 0};
    protected long[] nCheck = {0, 0};

    public ClusterWheatleySoft(int nPoints, MayerFunction f, double tol) {
        this.n = nPoints;
        this.f = f;
        int nf = 1 << n;  // 2^n
        fQ = new double[nf];
        fC = new double[nf];
        for (int i = 0; i < n; i++) {
            fQ[1 << i] = 1.0;
        }
        fA = new double[nf];
        fB = new double[nf];
        this.tol = tol;
        clusterBD = tol == 0 ? null : new ClusterWheatleySoftBD(nPoints, f, -3 * (int) Math.log10(tol));
    }

    public void setDoCaching(boolean newDoCaching) {
        doCaching = newDoCaching;
        if (clusterBD != null) {
            clusterBD.setDoCaching(doCaching);
        }
    }

    /**
     * Directs this cluster to only compute p fraction of the time when the
     * value is too small (below tol).  When it is computed, the value will be
     * boosted by 1/p.
     *
     * @param p   the fraction of time BD values will be computed
     * @param rng the random number generated used to decide to do BD or not
     */
    public void setBDAccFrac(double p, IRandom rng) {
        BDAccFrac = p;
        random = rng;
    }

    public ClusterAbstract makeCopy() {
        ClusterWheatleySoft c = new ClusterWheatleySoft(n, f, tol);
        c.setTemperature(1 / beta);
        c.setDoCaching(doCaching);
        if (BDAccFrac < 1) {
            c.setBDAccFrac(BDAccFrac, random);
        }
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
            debug = true;
            updateF(box);
            calcValue(box);
            throw new RuntimeException("oops "+value);
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
            if (fQ[i] == 0) continue;
            //loop over pairs formed from j and each point in partition; multiply by bond for each pair
            //all such pairs will be with bits higher than j, as j is the lowest bit in i
            for (int l=(j<<1); l<i; l=(l<<1)) {
                if ((l&i)==0) continue; //l is not in partition
                fQ[i] *= fQ[l | j];
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
                fC[i] -= fC[j] * fQ[jComp];//for fQ, flip the bits on j; use only those appearing in i
            }
        }

        // find fA1
        for (int i=2; i<nf; i+=2) {
            // all even sets don't contain 1
            //fA[i] = 0;
            fB[i] = fC[i];
        }
        fA[1] = 0;
        fB[1] = fC[1];
        for (int i=3; i<nf; i+=2) {
            // every set will contain 1
            fA[i] = 0;
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
                fA[i] += fB[j] * fC[jComp|1];
            }
            fB[i] -= fA[i];//remove from B graphs that contain articulation point at 0
        }

        for (int v=1; v<n; v++) {
            int vs1 = 1<<v;
            for (int i=vs1+1; i<nf; i++) {
                fA[i] = 0;
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
                        fA[i] += fB[j] * (fB[jComp|vs1] + fA[jComp|vs1]);
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
                        if (j == i) break;
                        fA[i] += fB[j] * (fB[jComp | vs1] + fA[jComp | vs1]);
                    }
                }

                fB[i] -= fA[i];//remove from B graphs that contain articulation point at v
            }
        }

        totCount++;
        double bfac = (1.0 - n) / SpecialFunctions.factorial(n);
        if (Math.abs(fB[nf - 1]) < tol * 100 && fB[nf - 1] != 0) {
            double r = BDAccFrac < 1 ? random.nextDouble() : 1;
            boolean justChecking = Math.abs(fB[nf - 1]) > tol;
            if (justChecking) {
                r /= 1000;
            }
            // integrand is too small for recursion to compute accurately.  we ought to do
            // BD, but it's expensive.  only do BD BDAccFrac of the time.  If we do it, then
            // boost the returned value by 1/BDAccFrac to account for the missed configurations
            boolean doBD = clusterBD != null && (BDAccFrac == 1 || r < BDAccFrac);
            if (doBD) {
                countBD++;
                timeBD -= System.nanoTime();
                double valueBD = clusterBD.value(box);
                timeBD += System.nanoTime();
                if (justChecking) {
                    int idx = (int) Math.log10(Math.abs(fB[nf - 1]) / tol);
                    if (idx >= nCheck.length) idx = nCheck.length - 1;
                    nCheck[idx]++;
                    avgAbsCheck[idx] += (Math.abs(fB[nf - 1]) - avgAbsCheck[idx]) / nCheck[idx];
                    avgAbsCheckBD[idx] += (Math.abs(valueBD) / bfac - avgAbsCheckBD[idx]) / nCheck[idx];
                    value = bfac * fB[nf - 1];
                } else {
                    value = valueBD / BDAccFrac;
                }
            } else if (justChecking) {
                // value is near tol, but we decided to use double
                value = bfac * fB[nf - 1];
            } else {
                value = 0;
            }
        } else {
            value = bfac * fB[nf - 1];
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
        for(int i=0; i<n-1; i++) {
            for(int j=i+1; j<n; j++) {
                double ff = f.f(aPairs.getAPair(i,j),cPairs.getr2(i,j), beta);
                if (false && Double.isNaN(ff)) {
                    f.f(aPairs.getAPair(i,j),cPairs.getr2(i,j), beta);
                    throw new RuntimeException("oops");
                }
//                if (Math.abs(ff) < 1e-14) ff = 0;
                fQ[(1<<i)|(1<<j)] = ff+1;
            }
        }
    }

    public void setTemperature(double temperature) {
        beta = 1 / temperature;
        if (clusterBD != null) {
            clusterBD.setTemperature(temperature);
        }
    }

    public double[] getAverageCheck() {
        return avgAbsCheck;
    }

    public double[] getAverageCheckBD() {
        return avgAbsCheckBD;
    }
}
