/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial.wheatley;

import etomica.graph.model.Graph;
import etomica.graph.model.impl.GraphImpl;
import etomica.graph.operations.GraphOp.GraphOpNull;
import etomica.graph.operations.MaxIsomorph;
import etomica.graph.operations.MaxIsomorph.MaxIsomorphParameters;
import etomica.math.SpecialFunctions;
import etomica.virial.AtomPairSet;
import etomica.virial.BoxCluster;
import etomica.virial.CoordinatePairSet;
import etomica.virial.MayerFunction;
import etomica.virial.cluster.ClusterAbstract;

import java.util.HashSet;

/**
 * This class calculates the sum of all biconnected clusters using Wheatley's
 * recursive formulation.
 * 
 * @author David Kofke and Andrew Schultz 
 */
public class ClusterWheatleyHS implements ClusterWheatley {

    protected final int n;
    protected final MayerFunction f;
    
    protected final double[] fQ, fC;
    protected final double[] fA, fB;
    protected long cPairID = -1, lastCPairID = -1;
    protected double value, lastValue;
    protected double beta;
    protected final byte[] outDegree;
    protected final int[] fullBondMask;
    protected final boolean[] cliqueSet;
    protected int cliqueCount, eCliqueCount;
    protected boolean precalcQ = true;
    protected final int[] cliqueList, eCliqueList;

    public ClusterWheatleyHS(int nPoints, MayerFunction f) {
        this.n = nPoints;
        this.f = f;
        int nf = 1<<n;  // 2^n
        fQ = new double[nf];
        fC = new double[nf];
        for(int i=0; i<n; i++) {
            fQ[1<<i] = 1.0;
        }
        fA = new double[nf];
        fB = new double[nf];
        outDegree = new byte[n];
        fullBondMask = new int[nf];
        cliqueSet = new boolean[nf];
        cliqueList = new int[nf];
        eCliqueList = new int[nf];
    }

    public ClusterAbstract makeCopy() {
        ClusterWheatleyHS c = new ClusterWheatleyHS(n, f);
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
      
      calcValue(box, true);
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
        eCliqueCount = 0;
        // generate all partitions and compute product of e-bonds for all pairs in partition
        for (int i=3; i<nf; i++) {
            int j = i & -i;//lowest bit in i
            if (i==j) continue; // 1-point set
            int k = i&~j; //strip j bit from i and set result to k
            int jj = k&-k;
            if (k == jj) continue; // 2-point set; these fQ's were filled when bonds were computed, so skip
            int kk = i&~jj; // strip jj bit from i
            // fQ is either 0 or 1
            fQ[i] = fQ[k] * fQ[kk] * fQ[j|jj];
            if (fQ[i] == 0) continue;
            eCliqueList[eCliqueCount] = i;
            eCliqueCount++;
        }
    }

    public boolean checkConfig(BoxCluster box) {
        if (false && total>0 && total%1000000==0) {
            System.out.println("screened: "+((double)screened)/total+"     zero-not-screened: "+((double)(total-screened-notzero))/(total-screened));
        }
        updateF(box);
        total++;
        screened++;
        edgeCount = 0;

        for (int i=0; i<n; i++) {
            outDegree[i] = 0;
        }
        int nf = 1<<n;
        for (int i=0; i<nf; i++) {
            fullBondMask[i] = 0;
        }

        for (int i=0; i<n-1; i++) {
            for (int j=i+1; j<n; j++) {
                int k = (1<<i)|(1<<j);
                boolean fBond = (fQ[k] == 0); 
                cliqueSet[k] = fBond;
                if (fBond) {
                    outDegree[i]++;
                    outDegree[j]++;
                    edgeCount++;
                    fullBondMask[1<<i] |= 1<<j;
                    fullBondMask[1<<j] |= 1<<i;
                }
            }
        }

        int maxEdges = n*(n-1)/2;
        if (edgeCount==maxEdges) {
            screened--;
            cliqueCount = nf-n-(n*(n-1)/2);
            calcFullFQ(box);
            return true;
        }
        if (edgeCount==maxEdges-1) return false;
        for (int i=0; i<n; i++) {
            if (outDegree[i] < 2) {
                return false;
            }
        }

        // the existence of a clique separator indicates that the value for
        // this configuration is zero.  Loop through all sets, considering
        // each as a clique separator.
        cliqueCount = 0;
        for (int i=1; i<nf-3; i++) {
            int j = i & -i;//lowest bit in i
            if (i==j) {
                // 1-point set.  check as an articulation point
                if (cliqueCheck(i, nf-1)) return false;
                continue;
            }
            int k = i&~j; //strip j bit from i and set result to k
            int jj = k & -k; // 2nd lowest bit
            if (k == jj) {
                // 2-point set.  cliqueSet[i] was set in the above loop
                // over pairs.
                if (cliqueSet[i] && cliqueCheck(i, nf-1)) return false;
                continue;
            }

            int kk = i&~jj; // strip jj bit from i

            cliqueSet[i] = cliqueSet[k] && cliqueSet[kk] && cliqueSet[j|jj];

            if (!cliqueSet[i]) continue;

            // i is a clique
            if (cliqueCheck(i, nf-1)) {
                return false;
            }
            cliqueList[cliqueCount] = i;
            cliqueCount++;
        }
        
        if (precalcQ) {
            calcFullFQ(box);
        }

        screened--;
        return true;
    }

    /**
     * Checks to see if the set |clique| is a clique separator.  The graph
     * connectivity is described by the fullBondMask array.
     */
    protected boolean cliqueCheck(int clique, int nfm1) {
        int cComp = nfm1^clique;
        int m = cComp & -cComp;
        if (m==cComp) {
            // there is only one point not in the clique
            return false;
        }
        int seen = m | clique;

        while (true) {
            int newSeen = seen;
            while (m != 0) {
                int lb = m & -m;
                newSeen |= fullBondMask[lb];
                m = m ^ lb;
            }
            if (newSeen == seen) {
                // we didn't pick up any new points
                return true;
            }
            if (newSeen == nfm1) {
                // we found all points
                return false;
            }
            m = newSeen - seen;
            seen = newSeen;
        }
    }

    /**
     * Returns the cluster value for the given configuration.  You must call
     * doCheck(BoxCluster) before calling this method.
     */
    public double calcValue(BoxCluster box) {
        calcValue(box, false);
        return value;
    }

    /*
     * Computation of sum of biconnected diagrams.
     */
    protected void calcValue(BoxCluster box, boolean doCheck) {
        if (doCheck && !checkConfig(box)) {
            value = 0;
            return;
        }
        if (!precalcQ) {
            calcFullFQ(box);
        }

        int nf = 1<<n;
        // The sums below run over subsets j of i that contain a fixed set of bits (jBits), with j != i.
        // We write j = jBits | s and enumerate the proper submasks s of rest = i ^ jBits in ascending
        // order via s = (s - rest) & rest.  This visits only valid subsets, in the same order as the
        // original increment-and-skip loops, so the floating-point sums are unchanged.

        //Compute the fC's
        for(int i=1; i<nf; i++) {
            fC[i] = fQ[i];
            int iLowBit = i & -i;
            int rest = i ^ iLowBit;
            for (int s=0; s!=rest; s=(s-rest)&rest) {
                // j = iLowBit|s, jComp = rest^s is the subset of i complementing j
                fC[i] -= fC[iLowBit|s] * fQ[rest^s];
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
            //jBits has 1 and next lowest bit in i
            int rest = i ^ jBits;
            for (int s=0; s!=rest; s=(s-rest)&rest) {//sum over partitions of i containing jBits
                // j = jBits|s, jComp = rest^s is the subset of i complementing j
                fA[i] += fB[jBits|s] * fC[(rest^s)|1];
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
                if (iLowBit == vs1) {
                    // v is the lowest bit; jBits is v and the next lowest bit
                    int ii = i ^ iLowBit;
                    jBits = vs1 | (ii & -ii);
                }
                else {
                    // jBits is the lowest bit and v
                    jBits = iLowBit | vs1;
                }
                int rest = i ^ jBits;
                for (int s=0; s!=rest; s=(s-rest)&rest) {//sum over partitions of i containing jBits
                    // j = jBits|s, jComp = rest^s is the subset of i complementing j
                    int jCompv = (rest^s) | vs1;
                    fA[i] += fB[jBits|s] * (fB[jCompv] + fA[jCompv]);
                }

                fB[i] -= fA[i];//remove from B graphs that contain articulation point at v
            }
        }

        value = (1-n)*fB[nf-1]/SpecialFunctions.factorial(n);
        if (value != 0) {
            notzero++;
            // disable check above and then enable this to see if non-zero
            // configurations would be screened
            if (false && !checkConfig(box)) {
                Graph g = new GraphImpl((byte)n);
                for (int i=0; i<n-1; i++) {
                    for (int j=i+1; j<n; j++) {
                        if ((fullBondMask[i] & (1<<j)) != 0) {
                            g.putEdge((byte)i, (byte)j);
                        }
                    }
                }
//                MaxIsomorph maxIso = new MaxIsomorph();
//                MaxIsomorphParameters mip = new MaxIsomorphParameters(new GraphOpNull(), MaxIsomorph.PROPERTY_ALL);
//                g = maxIso.apply(g, mip);
                String s = g.getStore().toNumberString();
                System.out.println("**** oops thought this was zero: "+s);
                checkConfig(box);
            }
        }
        else if (false) {
            // enable this to see what zero-value configurations are not being
            // screened
            Graph g = new GraphImpl((byte)n);
            for (int i=0; i<n-1; i++) {
                for (int j=i+1; j<n; j++) {
                    if ((fullBondMask[i] & (1<<j)) != 0) {
                        g.putEdge((byte)i, (byte)j);
                    }
                }
            }
            MaxIsomorph maxIso = new MaxIsomorph();
            MaxIsomorphParameters mip = new MaxIsomorphParameters(new GraphOpNull(), MaxIsomorph.PROPERTY_ALL);
            g = maxIso.apply(g, mip);
            String s = g.getStore().toNumberString();
            if (!zeroMaps.contains(s)) {
                System.out.println(s+" is zero");
                zeroMaps.add(s);
                for (String ss : zeroMaps) {
                    System.out.print(ss+",");
                }
                System.out.println();
            }
        }
    }
    
    /**
     * Returns edgeCount (number of overlaps) of configuration passed to
     * checkConfig
     */
    public int getEdgeCount() {
        return edgeCount;
    }
    
    public int getCliqueCount() {
        return cliqueCount;
    }
    
    public int getECliqueCount() {
        return eCliqueCount;
    }
    
    public int[] getFullBondMask() {
        return fullBondMask;
    }
    
    public int[] getCliques() {
        return cliqueList;
    }

    public int[] getECliques() {
        return eCliqueList;
    }

    /**
     * Returns outDegee (number of bonds for each point) of the configuration
     * passed to checkConfig
     */
    public byte[] getOutDegree() {
        return outDegree;
    }

    protected long total, notzero, screened;
    protected int edgeCount;
    HashSet<String> zeroMaps = new HashSet<String>();

    protected void updateF(BoxCluster box) {
        CoordinatePairSet cPairs = box.getCPairSet();
        AtomPairSet aPairs = box.getAPairSet();

        f.setBox(box);
        // recalculate all f values for all pairs
        for(int i=0; i<n-1; i++) {
            for(int j=i+1; j<n; j++) {
                fQ[(1<<i)|(1<<j)] = f.f(aPairs.getAPair(i,j),cPairs.getr2(i,j), beta)+1;
            }
        }
    }

    public void setTemperature(double temperature) {
        beta = 1/temperature;
    }

    public static class AllSigs extends ClusterWheatleyHS {
        protected long sig;
        
        public AllSigs(int nPoints) {
            super(nPoints, null);
        }
        
        public long valueForSig(long s) {
            sig = s;
            calcValue(null, true);
            return (long)value;
        }
        
        public void updateF(BoxCluster box) {
            long b = 1;
            for (int i=0; i<n-1; i++) {
                for (int j=i+1; j<n; j++) {
                    fQ[(1<<i)|(1<<j)] = ((sig & b) == 0) ? 1 : 0;
                    b <<= 1;
                }
            }
        }
    }
}
