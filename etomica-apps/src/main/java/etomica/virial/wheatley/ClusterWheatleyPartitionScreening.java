/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial.wheatley;

import etomica.graph.model.Graph;
import etomica.graph.model.impl.GraphImpl;
import etomica.graph.operations.GraphOp.GraphOpNull;
import etomica.graph.operations.MaxIsomorph;
import etomica.graph.operations.MaxIsomorph.MaxIsomorphParameters;
import etomica.virial.AtomPairSet;
import etomica.virial.BoxCluster;
import etomica.virial.CoordinatePairSet;
import etomica.virial.MayerFunction;
import etomica.virial.cluster.ClusterAbstract;

import java.util.Arrays;
import java.util.BitSet;
import java.util.HashSet;

/**
 * Computes the hard-sphere virial-coefficient integrand using Wheatley's recursion for
 * the sum over all biconnected Mayer diagrams, in exact integer arithmetic, with
 * screening of zero-valued configurations and lookup tables for small subsets.
 *
 * <h2>What is computed</h2>
 * For hard spheres every Mayer bond is either &minus;1 (overlap) or 0, so the sum
 * <i>f</i><sub>B</sub> over biconnected diagrams on <i>n</i> points is an integer for
 * each configuration. {@link #value(BoxCluster)} returns
 * <blockquote>
 *   (1 &minus; <i>n</i>) <i>f</i><sub>B</sub>.
 * </blockquote>
 * <b>Normalization:</b> this is <i>n</i>! times the value returned by
 * {@link ClusterWheatleyHS} for the same configuration, which returns
 * (1 &minus; <i>n</i>) <i>f</i><sub>B</sub>/<i>n</i>!. A driver that switches between the
 * two classes must account for the factor of <i>n</i>!.
 *
 * <h2>Algorithm</h2>
 * The recursion <i>f</i><sub>Q</sub> &rarr; <i>f</i><sub>C</sub> (connected) &rarr;
 * <i>f</i><sub>B</sub> (biconnected) follows R. J. Wheatley, Phys. Rev. Lett.
 * <b>110</b>, 200601 (2013). It differs from {@link ClusterWheatleyHS} in four ways:
 * <ul>
 * <li><b>Integer arithmetic.</b> {@code fC}, {@code fA} and {@code fB} are {@code int}
 *     arrays, so results are exact. The constructor rejects <i>n</i> &gt; 13, where
 *     {@code int} could overflow.</li>
 * <li><b>Boolean f<sub>Q</sub>.</b> {@code fQ[S]} is true when some pair in subset
 *     <i>S</i> overlaps, that is, when the product of e-bonds over <i>S</i> is 0.</li>
 * <li><b>Precomputed partitions.</b> The constructor lists, for every subset, the
 *     sub-subsets that each step of the recursion sums over, so each configuration
 *     visits only valid terms. This uses O(3<sup><i>n</i></sup>) memory.</li>
 * <li><b>Lookup tables.</b> For every subset of at most {@code nPtsTabulated} points,
 *     the constructor tabulates f<sub>C</sub>, f<sub>A</sub> and f<sub>A</sub> +
 *     f<sub>B</sub> (at each articulation-point level) for every possible pattern of
 *     overlaps (the subset's "signature"). For each configuration,
 *     {@link #calcSignatures()} finds each small subset's signature and its values are
 *     looked up instead of computed. If {@code nPtsTabulated} equals <i>n</i>, the
 *     result is read directly from the table. The tables have
 *     2<sup><i>p</i>(<i>p</i>&minus;1)/2</sup> entries for <i>p</i> =
 *     {@code nPtsTabulated}, so construction time and memory grow very rapidly with
 *     <i>p</i>.</li>
 * </ul>
 *
 * <h2>Screening</h2>
 * {@link #value(BoxCluster)} first calls {@link #checkConfig(BoxCluster)}, which builds
 * the overlap graph and returns false (value 0) without running the recursion if
 * <ul>
 * <li>exactly one pair does not overlap;</li>
 * <li>any point overlaps fewer than two others; or</li>
 * <li>the overlap graph has a clique separator: a clique (a single point, an
 *     overlapping pair, or a larger set in which every pair overlaps) whose removal
 *     disconnects the remaining points.</li>
 * </ul>
 * If every pair overlaps, the configuration is accepted without the clique search.
 * {@code checkConfig} also records the statistics returned by {@link #getEdgeCount()},
 * {@link #getOutDegree()}, {@link #getFullBondMask()}, {@link #getCliques()} and
 * {@link #getECliques()}.
 *
 * <h2>Inputs and conventions</h2>
 * <ul>
 * <li>{@code f} must return &minus;1 for an overlapping pair; any other value is
 *     treated as no overlap. {@link etomica.virial.MayerHardSphere} is the usual
 *     choice.</li>
 * <li>Values are cached by the configuration ID of the box's {@link CoordinatePairSet},
 *     including the configuration before the most recent one, so that a rejected Monte
 *     Carlo trial does not trigger a recomputation.</li>
 * </ul>
 *
 * <h2>Limitations</h2>
 * <ul>
 * <li>{@link #makeCopy()} uses the default {@code nPtsTabulated} of 2, not this
 *     cluster's value.</li>
 * <li>{@link #main(String[])}, {@link #updateF()}, the {@code checkme} and
 *     {@code doStatistics} flags and {@link FrequencyCounter} are development aids.
 *     Both flags are compile-time false, so {@link #sigCounter} is null.</li>
 * </ul>
 *
 * <p>Drivers: {@link etomica.virial.simulations.hardsphere.VirialHS} and
 * {@link etomica.virial.simulations.hardsphere.VirialHSBinMultiThreaded}, when run with
 * {@code nPtsTabulated} &gt; 0 (otherwise they use {@link ClusterWheatleyHS}).
 *
 * @author David Kofke and Andrew Schultz
 * @see ClusterWheatleyHS
 */
public class ClusterWheatleyPartitionScreening implements ClusterWheatley {

    /** Number of points (n) and number of subsets (2<sup>n</sup>). */
    protected final int n, nf;
    /** Hard-sphere Mayer function; &minus;1 means overlap. */
    protected final MayerFunction f;
    
    /** For each subset bitmask, true if some pair in the subset overlaps (product of e-bonds is 0). */
    protected final boolean[] fQ;
    /**
     * Recursion values for each subset bitmask: f<sub>A</sub> (articulation-point
     * diagrams at the current level), f<sub>B</sub>, f<sub>A</sub> + f<sub>B</sub> and
     * f<sub>C</sub> (connected diagrams).
     */
    protected final int[] fA, fB, fAB, fC;
    /** Configuration IDs for the current and previous cached values. */
    protected long cPairID = -1, lastCPairID = -1;
    /** Current and previous cached values. */
    protected double value, lastValue;
    /** Reciprocal temperature, passed to the Mayer function. */
    protected double beta;
    /** Number of overlaps at each point, from the most recent {@link #checkConfig(BoxCluster)}. */
    protected final byte[] outDegree;
    /**
     * For each single-point bitmask 1&lt;&lt;i, the bitmask of points that overlap point
     * i, from the most recent {@link #checkConfig(BoxCluster)}. Other entries are 0.
     */
    protected final int[] fullBondMask;
    /** For each subset bitmask, whether it is a clique of overlaps (work array for screening). */
    protected final boolean[] cliqueSet;
    /** Numbers of overlap cliques and of e-bond cliques from the most recent {@link #checkConfig(BoxCluster)}. */
    protected int cliqueCount, eCliqueCount;
    /** If true (always, at present), f<sub>Q</sub> is computed during {@link #checkConfig(BoxCluster)}. */
    protected boolean precalcQ = true;
    /** Lists of overlap cliques and e-bond cliques, as subset bitmasks. */
    protected final int[] cliqueList, eCliqueList;
    
//    protected final int[] nBonds;
    /** Number of points in each subset, indexed by subset bitmask. */
    protected final int[] nPts;
    /** Work array: articulation-point level reached for each tabulated subset. */
    protected final int[] vCount;
    /** Overlap signature of each tabulated subset: one bit per pair, set if the pair overlaps. */
    protected final int[] sig;
    /**
     * Precomputed partitions for the f<sub>A</sub> sums: {@code partitionsA[v][i]} lists,
     * in ascending order, the proper subsets of i that contain point v and the lowest
     * point of i other than v; null if i does not contain v or has fewer than three
     * points.
     */
    protected final int[][][] partitionsA;
    /**
     * Precomputed partitions for the f<sub>C</sub> sums: {@code partitionsC[i]} lists, in
     * ascending order, the proper subsets of i that contain the lowest point of i.
     */
    protected final int[][] partitionsC;
    long bigSum1, bigSum2, count, bigSumN, bigSumNCount;

    /** Tabulated f<sub>A</sub> values: {@code fAList[p][signature][v]} for subsets of p points. */
    protected final int[][][] fAList;
    /** Tabulated f<sub>A</sub> + f<sub>B</sub> values: {@code fABList[p][signature][v]} for subsets of p points. */
    protected final int[][][] fABList; 
    /** Tabulated f<sub>C</sub> values: {@code fCList[p][signature]} for subsets of p points. */
    protected final int[][] fCList;
    /** Shared all-zero table row, used for signatures with too few overlaps to be connected. */
    protected final int[] allZero;
    
    /** Table rows looked up for each subset of the current configuration; null if not tabulated. */
    protected final int[][] fAValues, fABValues;
    /** Subsets of at most this many points are tabulated (at most n). */
    protected final int nPtsTabulated;//tabulate fA,fAB,fC values for all graphs up to this size
    /** Development flag: if true, tabulated values are checked against computed ones. */
    protected static final boolean checkme = false;
    /** Development flag: if true, signature frequencies and maximum magnitudes are recorded. */
    protected final boolean doStatistics = false;
    /** Signature frequency counters, one per subset size; null unless {@code doStatistics} is true. */
    public final FrequencyCounter[] sigCounter;
    long maxA, maxB, maxAB, maxC, maxQ;

    /**
     * Creates the cluster, tabulating subsets of up to two points.
     *
     * @param nPoints number of points, n, at most 13
     * @param f       hard-sphere Mayer function (&minus;1 for overlap)
     * @throws IllegalArgumentException if nPoints &gt; 13
     */
    public ClusterWheatleyPartitionScreening(int nPoints, MayerFunction f) {
        this(nPoints, f, 2);
    }

    /**
     * Creates the cluster, precomputing the partition lists and the lookup tables.
     *
     * @param nPoints       number of points, n, at most 13
     * @param f             hard-sphere Mayer function (&minus;1 for overlap)
     * @param nPtsTabulated tabulate subsets of up to this many points (reduced to n if
     *                      larger); the tables have 2<sup>p(p&minus;1)/2</sup> entries
     *                      for p = nPtsTabulated
     * @throws IllegalArgumentException if nPoints &gt; 13
     */
    public ClusterWheatleyPartitionScreening(int nPoints, MayerFunction f, int nPtsTabulated) {
        this.n = nPoints;
        if(nPoints > 13) throw new IllegalArgumentException("fA, fB, etc need to be long instead of int for nPoints > 13");
        this.f = f;
        this.nPtsTabulated = nPoints > nPtsTabulated ? nPtsTabulated : nPoints;
        nf = 1<<n;  // 2^n
        fQ = new boolean[nf];
        fC = new int[nf];
        nPts = new int[nf];
        vCount = new int[nf];
        fAValues = new int[nf][];
        fABValues = new int[nf][];
        sig = new int[nf];
        allZero = new int[n];

        if(doStatistics) {
            sigCounter = new FrequencyCounter[n];
            for(int i=1; i<=n; i++) {
                sigCounter[i-1] = new FrequencyCounter(i, Math.min(2000, Math.abs(1<<nPairs(i))));
            }
        } else {
            sigCounter = null;
        }
        
        for(int i=1; i<nf; i++) {
            nPts[i] = Integer.bitCount(i);
        }
        for(int i=0; i<n; i++) {
            fQ[1<<i] = false;
        }
        fA = new int[nf];
        fB = new int[nf];
        fAB = new int[nf];
        outDegree = new byte[n];
        fullBondMask = new int[nf];
        cliqueSet = new boolean[nf];
        cliqueList = new int[nf];
        eCliqueList = new int[nf];

        partitionsC = new int[nf][];
        for(int i=1; i<nf; i++) {
            partitionsC[i] = computeCPartitions(i);
        }

        partitionsA = new int[n][nf][];
        computeAPartitions();
        
        fAList = new int[n+1][][];
        fABList = new int[n+1][][];
        fCList = new int[n+1][];
        computefTables();
        
        if(checkme) System.out.println("Note -- checkme is true");
//        System.out.println("tabulation complete");
//        System.exit(1);
        
        if(doStatistics) {
            System.out.println("Note -- doStatistics is true");
            for(int i=0; i<n; i++) {
                sigCounter[i].reset();
            }
        }        
   }
    
    private final int[] computeCPartitions(int i) {
        // proper subsets j of i that contain the lowest bit of i, in ascending order
        int iLowBit = i & -i;
        int rest = i ^ iLowBit;
        int[] iPartitions = new int[(1 << Integer.bitCount(rest)) - 1];
        int k = 0;
        for (int s = 0; s != rest; s = (s - rest) & rest) {
            iPartitions[k++] = iLowBit | s;
        }
        return iPartitions;
    }
    
    private final void computeAPartitions() {
        for (int v=0; v<n; v++) {
            int vs1 = 1<<v;
            int[][] vPartitions = partitionsA[v];
            for (int i=vs1+1; i<nf; i++) {
                vPartitions[i] = computeAPartitions(vs1, i);
            }
        }
    }
        
    /**
     * Lists the partitions summed over when computing f<sub>A</sub> for subset i at
     * articulation point v.
     *
     * @param vs1 the bitmask of point v (1&lt;&lt;v)
     * @param i   the subset bitmask
     * @return the proper subsets of i that contain v and the lowest point of i other
     *         than v, in ascending order; null if i does not contain v or has fewer than
     *         three points
     */
    protected int[] computeAPartitions(int vs1, int i) {
        
        if((vs1 & i)==0 || (nPts[i]<3)) return null;

        // proper subsets j of i that contain v and the lowest bit of i other than v,
        // in ascending order
        int iLowBit = (i&-i);//lowest bit in i
        int jBits;
        if (iLowBit == vs1) {
            int ii = i ^ iLowBit;
            jBits = vs1 | (ii & -ii);
        }
        else {
            jBits = iLowBit | vs1;
        }
        int rest = i ^ jBits;
        int[] iPartitions = new int[(1 << Integer.bitCount(rest)) - 1];
        int k = 0;
        for (int s = 0; s != rest; s = (s - rest) & rest) {
            iPartitions[k++] = jBits | s;
        }
        return iPartitions;
    }
    
    
    

    private final void computefTables() {
        for(int np=2; np<=nPtsTabulated; np++) {
            fAList[np] = new int[1<<nPairs(np)][np];
            fABList[np] = new int[1<<nPairs(np)][np];
            fCList[np] = new int[1<<nPairs(np)];
        }
        int sigMax = (1 << nPairs(nPtsTabulated));
        long nDuplicate = 0;
        long nArrays = 0;
//        long t1 = System.currentTimeMillis();
        //construct table of fA, fB, and fC for each signature for each nPts being tabulated
        for(int s=0; s<sigMax; s++) {//loop over signatures
            int nb = Integer.bitCount(s);
//            if (s>0 && s%10000==0) {
//                long t2 = System.currentTimeMillis();
//                System.out.println(String.format("s=%d/%d  time: %d/%d", s, sigMax, (t2-t1)/1000, (int)(((double)(t2-t1))*sigMax/s/1000)));
//            }
            for(int np=2; np<=nPtsTabulated; np++) {
                if(s >= (1<<nPairs(np))) continue;//signature has more bonds than can be formed by np points
                if(nb+1 < np) {
                    nDuplicate += 2;
                    nArrays++;
                    fAList[np][s] = allZero;
                    fABList[np][s] = allZero;
                    fCList[np][s] = 0;
                    continue;
                }
                //construct bonds from signature
                int bondBit = 1;
                for(int j=(1<<(np-2)); j>0; j=(j>>1)) {
                    for (int l=(j<<1); l<(1<<np); l=(l<<1)) {
                        fQ[l | j] = ((s & bondBit) == bondBit);
                        bondBit = bondBit << 1;
                    }
                }

                calcfQ(null);

//                System.out.println(np+"\t"+s+"\t"+sig[(1<<np)-1]+"\t"+(s-sig[(1<<np)-1]));
            
                calcfC(false);
                int i = (1<<np) - 1;
                fCList[np][s] = fC[i];
                for(int v=0; v<np; v++) {
                    calcfAB(v, false);
                    fAList[np][s][v] = fA[i];
                    fABList[np][s][v] = fAB[i];
                }//end of v-loop

                //remove duplicate arrays to save memory
                nArrays++;
                for(int s0=0; s0<s; s0++) {
                    if(Arrays.equals(fAList[np][s0], fAList[np][s])) {
                        fAList[np][s] = fAList[np][s0];
                        nDuplicate++;
                        break;
                    }
                }
                for(int s0=0; s0<s; s0++) {
                    if(Arrays.equals(fABList[np][s0], fABList[np][s])) {
                        fABList[np][s] = fABList[np][s0];
                        nDuplicate++;
                        break;
                    }
                }
 
            }//end of np-loop
        }//end of s-loop
        
//        System.out.println("Removed "+nDuplicate+" duplicate of "+(2*nArrays)+" total arrays");
    }

    private static int nPairs(int nPts) {
        return nPts*(nPts-1)/2;
    }

    /**
     * Returns a new cluster with the same points, Mayer function and temperature. Note
     * that the copy tabulates subsets of up to two points, whatever this cluster uses.
     */
    public ClusterAbstract makeCopy() {
        ClusterWheatleyPartitionScreening c = new ClusterWheatleyPartitionScreening(n, f);
        c.setTemperature(1/beta);
        return c;
    }

    /** Returns the number of points, n. */
    public int pointCount() {
        return n;
    }

    /**
     * Returns (1 &minus; n) f<sub>B</sub> for the box's current configuration, or 0 if it
     * is screened out. See the class description for the normalization, which differs
     * from {@link ClusterWheatleyHS} by a factor of n!.
     *
     * @param box the box holding the configuration
     * @return the integer-valued cluster value
     */
    public double value(BoxCluster box) {
      CoordinatePairSet cPairs = box.getCPairSet();
      long thisCPairID = cPairs.getID();
      if (thisCPairID == cPairID) {
          return value;
      }
      if(thisCPairID == lastCPairID) {
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
      
      calcValue(box, true);
      if (Double.isNaN(value) || Double.isInfinite(value)) {
          updateF(box);
          calcValue(box);
          throw new RuntimeException("oops");
      }
      return value;
    }

    /** {@return the configuration ID of the current cached value} */
    public long getCPairID() {
        return cPairID;
    }
    
    /** {@return the configuration ID of the previous cached value} */
    public long getLastCPairID() {
        return lastCPairID;
    }

    /** {@return the previous cached value} */
    public double getLastValue() {
        return lastValue;
    }

    /**
     * Fills {@link #fQ} for every subset of three or more points, given that the pair
     * entries have already been set by {@link #updateF(BoxCluster)}, and lists the
     * e-bond cliques (subsets with no overlapping pair).
     *
     * @param box the box holding the configuration (not used)
     */
    protected void calcfQ(BoxCluster box) {
        eCliqueCount = 0;
        // generate all partitions and compute product of e-bonds for all pairs in partition
loop1:  for (int i=7; i<nf; i++) {
            if(nPts[i] < 3) continue;//1-point sets filled on construction; 2-point set fQ's were filled when bonds were computed
            int j = i & -i;//lowest bit in i
            int k = i&~j; //strip j bit from i and set result to k
            fQ[i] = fQ[k]; //initialize with previously-computed product of all pairs in partition, other than j
            if(fQ[i]) continue;//nothing will change if fQ[i] is or becomes true

            //loop over pairs formed from j and each point in partition; multiply by bond for each pair
            //all such pairs will be with bits higher than j, as j is the lowest bit in i
            for (int l=(j<<1); l<i; l=(l<<1)) {
                if ((l&i)==0) continue; //l is not in partition
                //fQ[i] |= fQ[l|j], but can make equality here because loop breaks if fQ[i] is true
                fQ[i] = fQ[l | j];//fQ = true indicates a non-zero f-bond (e = 0); so product of e-bonds will be 0 (fQ[i] = true) if any of them are zero
                if(fQ[i]) continue loop1;
            }
            eCliqueList[eCliqueCount] = i;
            //if we get here, fQ[i] is false
            eCliqueCount++;//for fQ[i]=false, all e-bonds are non-zero
        }
    }
    
    /**
     * Computes the overlap signature of every subset of at most {@code nPtsTabulated}
     * points and looks up its tabulated f<sub>A</sub>, f<sub>A</sub> + f<sub>B</sub> and
     * f<sub>C</sub> values.
     */
    protected void calcSignatures() {
        
        for (int i=3; i<nf; i++) {
            if(((nPts[i] > nPtsTabulated) || (nPts[i] == 1)) && !doStatistics) continue;

            if(nPts[i] == 2) {
                sig[i] = fQ[i] ? 1 : 0;
                fAValues[i] = fAList[2][sig[i]];
                fABValues[i] = fABList[2][sig[i]];
                fC[i] = fCList[2][sig[i]];
                if(doStatistics) sigCounter[1].add(sig[i]);
                continue;               
            }
            
            int j = i & -i;//lowest bit in i
            int k = i^j; //strip j bit from i and set result to k
            sig[i] = sig[k];
            int bondBit = 1 << nPairs(nPts[k]);
            //loop over bits in i, excluding j
            for(int it=k, l=(it&-it); it>0; it^=l,l=(it&-it)) {
                if(fQ[l | j]) {
                    sig[i] |= bondBit;
                }
                bondBit = bondBit << 1;
            }
            
            fAValues[i] = fAList[nPts[i]][sig[i]];
            fABValues[i] = fABList[nPts[i]][sig[i]];
            fC[i] = fCList[nPts[i]][sig[i]];
            
            if(doStatistics) {
                if(maxC < Math.abs(fC[i])) {
                    maxC = Math.max(maxC, Math.abs(fC[i]));
                    System.out.println("MaxA,B,AB,C: "+maxA+" "+maxB+" "+maxAB+" "+maxC);
                }
                sigCounter[nPts[i]-1].add(sig[i]);
            }
        }
    }


    /**
     * Evaluates the Mayer function (via {@link #updateF(BoxCluster)}), builds the overlap
     * graph and applies the screening tests described in the class description. If the
     * configuration is not screened out, also computes f<sub>Q</sub> and records the
     * statistics returned by the getters.
     *
     * @param box the box holding the configuration
     * @return false if the value is zero for this configuration; true otherwise
     */
    public boolean checkConfig(BoxCluster box) {
        if (false && total>0 && total%100000==0) {
            System.out.println("screened: "+((double)screened)/total+"     zero-not-screened: "+((double)(total-screened-notzero))/(total-screened));
        }
        updateF(box);
        total++;
        screened++;
        edgeCount = 0;

        for (int i=0; i<n; i++) {
            outDegree[i] = 0;
        }

        for (int i=0; i<nf; i++) {
            fullBondMask[i] = 0;
        }

        for (int i=0; i<n-1; i++) {
            for (int j=i+1; j<n; j++) {
                int k = (1<<i)|(1<<j);
                boolean fBond = fQ[k]; 
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
            calcfQ(box);
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
iLoop:  for (int i=1; i<nf-3; i++) {
            int j = i & -i;//lowest bit in i
            if (i==j) {
                // 1-point set.  check as an articulation point
                if (cliqueCheck(i, nf-1)) return false;
                continue;
            }
            int k = i&~j; //strip j bit from i and set result to k
            if (k == (k&-k)) {
                // 2-point set.  cliqueSet[i] was set in the above loop
                // over pairs.
                if (cliqueSet[i] && cliqueCheck(i, nf-1)) return false;
                continue;
            }
            cliqueSet[i] = cliqueSet[k]; //initialize with previously-computed product of all pairs in partition, other than j

            if (!cliqueSet[i]) continue;
            //loop over pairs formed from j and each point in partition; multiply by bond for each pair
            //all such pairs will be with bits higher than j, as j is the lowest bit in i
            for (int l=(j<<1); l<i; l=(l<<1)) {
                if ((l&i)==0) continue; //l is not in partition
                if (!cliqueSet[l|j]) {
                    cliqueSet[i] = false;
                    continue iLoop;
                }
            }
            // i is a clique
            if (cliqueCheck(i, nf-1)) {
                return false;
            }
            cliqueList[cliqueCount] = i;
            cliqueCount++;
        }
        
        if (precalcQ) {
            calcfQ(box);
        }

        screened--;
        return true;
    }

    /**
     * Checks whether a set is a separator: whether removing it leaves the remaining
     * points disconnected in the overlap graph described by {@link #fullBondMask}.
     *
     * @param clique the bitmask of the set to test
     * @param nfm1   the bitmask of all points, 2<sup>n</sup> &minus; 1
     * @return true if the set separates the graph
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
     * Runs the recursion and returns the value, without screening.
     * {@link #checkConfig(BoxCluster)} must have been called first for this
     * configuration (it computes f<sub>Q</sub>). Most callers should use
     * {@link #value(BoxCluster)}.
     *
     * @param box the box holding the configuration
     * @return (1 &minus; n) f<sub>B</sub>
     */
    public double calcValue(BoxCluster box) {
        calcValue(box, false);
        return value;
    }

    /**
     * Computes the value and stores it in {@link #value}.
     *
     * @param box     the box holding the configuration
     * @param doCheck if true, call {@link #checkConfig(BoxCluster)} first and set the
     *                value to 0 if the configuration is screened out
     */
    protected void calcValue(BoxCluster box, boolean doCheck) {
        if (doCheck && !checkConfig(box)) {
            value = 0;
            return;
        }
        if (!precalcQ) {
            calcfQ(box);
        }
        
        calcSignatures();
        
        if(fAValues[nf-1] != null) {//desired value is tabulated, so bypass whole calculation
            value = (1-n)*(fABValues[nf-1][n-1] - fAValues[nf-1][n-1]);
            
        } else {

            calcfC(true);
            
            for(int i=1; i<nf; i++) {
                vCount[i] = 0;
            }
    
            for (int v=0; v<n; v++) {
                calcfAB(v, true);
            }
        
            value = (1-n)*fB[nf-1];
        }

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
     * Computes f<sub>C</sub> for every subset.
     *
     * @param useTable if true, skip subsets whose value was looked up by
     *                 {@link #calcSignatures()}
     */
    protected final void calcfC(final boolean useTable) {

        for(int i=1; i<nf; i++) {
            if((fAValues[i] != null) && useTable) {
                continue;//fC was set in calcSignatures
            }
            
            fC[i] = fQ[i] ? 0 : 1;

            int[] iPartitions = partitionsC[i];
            final int kmax = iPartitions.length;
            for(int k=0; k<kmax; k++) {
                int j = iPartitions[k];
                int jComp = (i & ~j);
                if(!fQ[jComp]) fC[i] -= fC[j];
            }//end k-loop
                        
            if(doStatistics) {
                if(maxC < Math.abs(fC[i])) {
                    maxC = Math.max(maxC, Math.abs(fC[i]));
                    System.out.println("MaxA,B,AB,C: "+maxA+" "+maxB+" "+maxAB+" "+maxC);
                }
            }
        }
    }
    
    /**
     * Applies articulation-point level v of the recursion: computes f<sub>A</sub>
     * (diagrams with an articulation point at v but none below) and removes it from
     * f<sub>B</sub>, for every subset.
     *
     * @param v        the articulation point, 0 &le; v &lt; n
     * @param useTable if true, use tabulated values where available
     */
    protected final void calcfAB(int v, final boolean useTable) {
        
        if(v == 0) {
            if(checkme && useTable) {
                calcfAB0(false);
                for(int i=1; i<nf; i++) {
                    if(fAValues[i] != null) {
                        int fATable = (i%2==1) ? fAValues[i][vCount[i]] : 0;
                        int fABTable = (i%2==1) ? fABValues[i][vCount[i]++] : fC[i];
                        int fBTable = fABTable - fATable;
                        int fAerr = fA[i]-fATable;
                        int fBerr = fB[i]-fBTable;
                        int fABerr = fAB[i]-fABTable;
                        if(fAerr != 0 || fABerr != 0) System.out.println(v+"\t"+i+"\t"+ fA[i] +"\t"+fATable+"\t"+fAerr+"\t"+ fAB[i] +"\t"+fABTable+"\t"+fABerr+"\t"+ fB[i] +"\t"+fBTable+"\t"+fBerr);
                    }
                }
            } else {
                calcfAB0(useTable);
            }
            return;
        }//end v==0

        int vs1 = 1<<v;
        int[][] vPartitions = partitionsA[v];
        for (int i=vs1+1; i<nf; i++) {
            
            int[] iPartitions = vPartitions[i];            
            if (iPartitions == null) {
                //i doesn't contain v, or has too few vertices
                fA[i] = 0;
                fAB[i] = fB[i];
                //fB is unchanged
                continue;
                //need not concern about bypassing vCount increment, because
                //iPartitions==null means v is not in i, or nPts < 3
            }
            
            //fB[v][i] = fB[v-1][i];//no a.p. at v or below, starts with those having no a.p. at v-1 or below
            //rest of this is to generate A (diagrams having a.p. at v but not below), and subtract it from B
            if(useTable && (fAValues[i] != null) && !checkme) {
                fA[i] = fAValues[i][vCount[i]];
                fAB[i] = fABValues[i][vCount[i]];
                vCount[i]++;
            } else {
                fAB[i] = fB[i];//fAB(v) = fB(v-1) because fB(v) = fB(v-1) - fA(v)
                fA[i] = 0;//initialize fA(v)
                final int kmax = iPartitions.length;
                for(int k=0; k<kmax; k++) {
                    int j = iPartitions[k];
                    if(fB[j]==0) continue;//this is true about 90% of the time
                    int jComp = (i & ~j) | vs1;
                    fA[i] += fB[j] * fAB[jComp];
                }//end k-loop

                //this code is used just when checking that tables are working correctly
                if(checkme && (fAValues[i] != null) & useTable) {
                    int fATable = fAValues[i][vCount[i]];
                    int fABTable = fABValues[i][vCount[i]]; 
                    int fBTable = fABTable - fATable;
                    int fAerr = fA[i]-fATable;
                    int fBerr = fB[i]-fBTable;
                    int fABerr = fAB[i]-fABTable;
                    if(fAerr != 0 || fABerr != 0) System.out.println(v+"\t"+i+"\t"+ fA[i] +"\t"+fATable+"\t"+fAerr+"\t"+ fAB[i] +"\t"+fABTable+"\t"+fABerr+"\t"+ fB[i] +"\t"+fBTable+"\t"+fBerr);
                    vCount[i]++;
                }//end if(checkme)
            }//end if
            fB[i] -= fA[i];//remove from B graphs that contain articulation point at v
            
            if(doStatistics) {
                if(maxA < Math.abs(fA[i]) || maxB < Math.abs(fB[i]) || maxAB < Math.abs(fAB[i])) {
                    maxA = Math.max(maxA, Math.abs(fA[i]));
                    maxB = Math.max(maxB, Math.abs(fB[i]));
                    maxAB = Math.max(maxAB,Math.abs(fA[i]+fB[i]));
                    System.out.println("MaxA,B,AB,C: "+maxA+" "+maxB+" "+maxAB+" "+maxC);
                }
            }
        }//end of i-loop
    }

    /**
     * Applies articulation-point level 0 of the recursion; see {@link #calcfAB(int, boolean)}.
     *
     * @param useTable if true, use tabulated values where available
     */
    protected final void calcfAB0(final boolean useTable) {

        for (int i=2; i<nf; i+=2) {// all even sets don't contain 1
            fA[i] = 0;
            fB[i] = fC[i];
            fAB[i] = fC[i];
        }

        fA[1] = 0;
        fB[1] = fC[1];
        fAB[1] = fC[1];
        int[][] vPartitions = partitionsA[0];
        for (int i=3; i<nf; i+=2) {// every set will contain 1
            if(fAValues[i] != null && useTable) {
                fA[i] = fAValues[i][0];
                vCount[i] = 1;
            } else {
                fA[i] = 0;
                int[] iPartitions = vPartitions[i];
                if(iPartitions != null) {
                    final int kmax = iPartitions.length;
                    for(int k=0; k<kmax; k++) {
                        int j = iPartitions[k];
                        int jComp = (i & ~j) | 1;
                        fA[i] += fB[j] * fC[jComp];
                    }
                }
            }
            fAB[i] = fC[i];
            fB[i] = fAB[i] - fA[i];//B is connected graphs that do not contain articulation point at 0

            if(doStatistics) {
                if(maxA < Math.abs(fA[i]) || maxB < Math.abs(fB[i]) || maxAB < Math.abs(fAB[i])) {
                    maxA = Math.max(maxA, Math.abs(fA[i]));
                    maxB = Math.max(maxB, Math.abs(fB[i]));
                    maxAB = Math.max(maxAB,Math.abs(fA[i]+fB[i]));
                    System.out.println("MaxA,B,AB,C: "+maxA+" "+maxB+" "+maxAB+" "+maxC);
                }
            }
        }

    }
    
    /** Returns the number of overlapping pairs found by the most recent {@link #checkConfig(BoxCluster)}. */
    public int getEdgeCount() {
        return edgeCount;
    }
    
    /**
     * Returns the number of overlap cliques found by the most recent
     * {@link #checkConfig(BoxCluster)}. When every pair overlaps, this is set to
     * 2<sup>n</sup> &minus; n &minus; n(n&minus;1)/2 and the list returned by
     * {@link #getCliques()} is not filled.
     */
    public int getCliqueCount() {
        return cliqueCount;
    }
    
    /**
     * Returns the number of e-bond cliques (subsets of three or more points with no
     * overlapping pair) found by the most recent {@link #checkConfig(BoxCluster)}.
     */
    public int getECliqueCount() {
        return eCliqueCount;
    }
    
    /** Returns {@link #fullBondMask}, an internal array. */
    public int[] getFullBondMask() {
        return fullBondMask;
    }
    
    /**
     * Returns the overlap cliques of three or more points as subset bitmasks. Only the
     * first {@link #getCliqueCount()} entries are valid. This is an internal array.
     */
    public int[] getCliques() {
        return cliqueList;
    }
    
    /**
     * Returns the e-bond cliques as subset bitmasks. Only the first
     * {@link #getECliqueCount()} entries are valid. This is an internal array.
     */
    public int[] getECliques() {
        return eCliqueList;
    }

    /**
     * Returns the number of overlaps at each point, from the most recent
     * {@link #checkConfig(BoxCluster)}. This is an internal array.
     */
    public byte[] getOutDegree() {
        return outDegree;
    }

    /** Counts of configurations checked, found non-zero, and screened out. */
    protected long total, notzero, screened;
    /** Number of overlapping pairs from the most recent {@link #checkConfig(BoxCluster)}. */
    protected int edgeCount;
    HashSet<String> zeroMaps = new HashSet<String>();

    /**
     * Evaluates the Mayer function for every pair and sets the pair entries of
     * {@link #fQ} (true for overlap).
     *
     * @param box the box holding the configuration
     */
    protected void updateF(BoxCluster box) {
        CoordinatePairSet cPairs = box.getCPairSet();
        AtomPairSet aPairs = box.getAPairSet();

        f.setBox(box);
        // recalculate all f values for all pairs
        for(int i=0; i<n-1; i++) {
            for(int j=i+1; j<n; j++) {
                fQ[(1<<i)|(1<<j)] = (f.f(aPairs.getAPair(i,j),cPairs.getr2(i,j), beta) == -1);
            }
        }
    }
    
    /**
     * Returns a bitmask with one bit per pair (in the order (0,1), (0,2), &hellip;,
     * (1,2), &hellip;), set if the Mayer function is non-zero for that pair. Used only
     * for debugging.
     *
     * @param box the box holding the configuration
     * @return the pair bitmask
     */
    public long calcSignature(BoxCluster box) {
        CoordinatePairSet cPairs = box.getCPairSet();
        AtomPairSet aPairs = box.getAPairSet();

        f.setBox(box);
        // recalculate all f values for all pairs
        int s = 0;
        long sig = 0;
        for(int i=0; i<n-1; i++) {
            for(int j=i+1; j<n; j++) {
                if (f.f(aPairs.getAPair(i,j),cPairs.getr2(i,j), beta) != 0) {
                  sig |= (1L<<s);
                }
                s++;
            }
        }
        return sig;
    }
    
    /** Development aid: sets a hand-chosen bond pattern instead of evaluating the Mayer function. */
    protected void updateF() {
        for(int i=0; i<n-1; i++) {
            for(int j=i+1; j<n; j++) {
                fQ[(1<<i)|(1<<j)] = true;//f.f(aPairs.getAPair(i,j),cPairs.getr2(i,j), beta)+1;
            }
        }
        fQ[(1<<2)|(1<<0)] = false;
        fQ[(1<<0)|(1<<1)] = false;
        fQ[(1<<1)|(1<<2)] = false;
    }

    /**
     * Sets the temperature passed to the Mayer function (hard-sphere results do not
     * depend on it).
     *
     * @param temperature the temperature, in simulation energy units (kT)
     */
    public void setTemperature(double temperature) {
        beta = 1/temperature;
    }
    
    /**
     * Development aid that records how often each overlap signature occurs for subsets of
     * a given size, to judge how much of the work tabulation would save. Used only when
     * {@code doStatistics} is true.
     */
    public final class FrequencyCounter {
        final int n;
        final int nfreq;
        int ncull = 100;
        final int[] sigs;
        final int[] freq;
        int counter = 0;
        int entries = 0;
        FrequencyCounter(int n, int nfreq) {
            if(nfreq < 0) nfreq = 0;
            if (ncull > nfreq/2) ncull=nfreq/2;
            this.n = n;
            this.nfreq = nfreq;
            sigs = new int[nfreq];
            freq = new int[nfreq];
            reset();
        }
        void reset() {
            for(int i=0; i<nfreq; i++) {
                sigs[i] = -1;
                freq[i] = 0;
                entries = 0;
                counter = 0;
            }
        }
        void add(int sig) {
            counter++;
            //look for instance
            for(int i=0; i<nfreq; i++) {
                if(sigs[i] == sig) {
                    freq[i]++;
                    return;
                }
            }
            //start a new instance
            if (entries == nfreq) {
                for(int k=0; k<ncull; k++) {
                    cull();
                }
            }
            for(int i=0; i<nfreq; i++) {
                if(sigs[i] < 0) {
                    sigs[i] = sig;
                    freq[i] = 1;
                    entries++;
                    return;
                }
            }
        }
        //method to remove low-frequency entries, if frequency array is filled
        void cull() {
            if(entries < nfreq-ncull) return;//no need to cull if there's still more empty slots than targeted buffer size
            int min = Integer.MAX_VALUE;
            int imin = 0;
            for(int i=0; i<nfreq; i++) {
                if(freq[i] < min && sigs[i] >= 0) {
                    imin = i;
                    min = freq[i];
                }
            }
            sigs[imin] = -1;
            freq[imin] = 0;
            entries--;
        }
        double[] getFractions() {
            double[] fracs = new double[nfreq];
            for(int i=0; i<nfreq; i++) {
                fracs[i] = (double)freq[i]/counter;
            }
            return fracs;
        }
        int[] getSigs() {
            return sigs;
        }
        /** Prints each signature that occurred more than once, with its frequency. */
        public void print() {
            System.out.println("Signature frequencies for n = "+n);
            double[] fractions = getFractions();
            for(int i=0; i<nfreq; i++) {
                if(freq[i]>1) System.out.println(String.format("%10d %7.5f %d", sigs[i], fractions[i], freq[i]));
            }
            System.out.println();
        }
        /** {@return the fraction of recorded occurrences held in the counter's table} */
        public double fractionTabulated() {
            int sum = 0;
            for(int i=0; i<nfreq; i++) {
                sum += freq[i];
            }
            return (double)sum/counter;
        }
        /** {@return the subset size this counter records} */
        public int getn() {
            return n;
        }
        /** {@return the number of distinct signatures held} */
        public int getEntries() {
            return entries;
        }
    }
    
    /**
     * Development aid: evaluates a hand-chosen bond pattern for six points.
     *
     * @param args not used
     */
    public static void main(String[] args) {
        int n = 6;
        BitSet bits = new BitSet(20);
        bits.set(5);
        bits.set(15);
        long i = (1<<5) | (1L<<63);
        System.out.println(Long.toBinaryString(i));
        System.out.println(bits.toString());
        ClusterWheatleyPartitionScreening cw = new ClusterWheatleyPartitionScreening(n, null);
        cw.updateF();
        cw.calcfQ(null);
        cw.calcValue(null, false);
//        for(int k=0; k<=n; k++) {
//            System.out.print(fAVals[k]);
//            if(k<n) System.out.print(", ");
//        }
//        System.out.println();
//        for(int k=0; k<=n; k++) {
//            System.out.print(fBVals[k]);
//            if(k<n) System.out.print(", ");
//        }
        
    }
}

