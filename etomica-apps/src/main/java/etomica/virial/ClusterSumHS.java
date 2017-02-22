/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial;


/**
 * This class calculates the sum of all biconnected clusters for hard spheres.
 * Configurations are screened, but the value is computed using standard
 * diagrams.
 * 
 * @author Andrew Schultz 
 */
public class ClusterSumHS extends ClusterSum implements ClusterWheatley {

    protected final int n;
    protected final double[] fQ, fC;
    protected final double[] fA, fB;
    protected final byte[] outDegree;
    protected final int[] fullBondMask;
    protected final boolean[] cliqueSet;
    protected int cliqueCount, eCliqueCount;
    protected boolean precalcQ = false;
    protected final int[] cliqueList, eCliqueList;

    public ClusterSumHS(ClusterBonds[] subClusters, double[] subClusterWeights, MayerFunction f) {
        super(subClusters, subClusterWeights, new MayerFunction[]{f});
        this.n = subClusters[0].pointCount();
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
        ClusterSumHS c = new ClusterSumHS(clusters, clusterWeights, f[0]);
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
        updateF(box);
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
            cliqueCount = nf-n-(n*(n-1)/2);
            if (precalcQ) calcFullFQ(box);
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
        return super.value(box);
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

    protected int edgeCount;

    protected void updateF(BoxCluster box) {
        CoordinatePairSet cPairs = box.getCPairSet();
        AtomPairSet aPairs = box.getAPairSet();

        f[0].setBox(box);
        // recalculate all f values for all pairs
        for(int i=0; i<n-1; i++) {
            for(int j=i+1; j<n; j++) {
                double fv = f[0].f(aPairs.getAPair(i,j),cPairs.getr2(i,j), beta);
                fQ[(1<<i)|(1<<j)] = fv + 1;
                fValues[i][j][0] = fv;
                fValues[i][j][1] = fv + 1;
            }
        }
    }

    public void setTemperature(double temperature) {
        beta = 1/temperature;
    }

}
