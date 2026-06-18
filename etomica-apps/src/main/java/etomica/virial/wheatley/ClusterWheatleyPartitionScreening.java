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
 * This class calculates the sum of all biconnected clusters using Wheatley's
 * recursive formulation.
 * 
 * @author David Kofke and Andrew Schultz 
 */
public class ClusterWheatleyPartitionScreening implements ClusterWheatley {

    protected final int n, nf;
    protected final MayerFunction f;
    
    protected final boolean[] fQ;
    protected final int[] fA, fB, fAB, fC;
    protected long cPairID = -1, lastCPairID = -1;
    protected double value, lastValue;
    protected double beta;
    protected final byte[] outDegree;
    protected final int[] fullBondMask;
    protected final boolean[] cliqueSet;
    protected int cliqueCount, eCliqueCount;
    protected boolean precalcQ = true;
    protected final int[] cliqueList, eCliqueList;
    
//    protected final int[] nBonds;
    protected final int[] nPts;
    protected final int[] vCount;
    protected final int[] sig;
    protected final int[][][] partitionsA;
    protected final int[][] partitionsC;
    long bigSum1, bigSum2, count, bigSumN, bigSumNCount;

    protected final int[][][] fAList;
    protected final int[][][] fABList; 
    protected final int[][] fCList;
    protected final int[] allZero;
    
    protected final int[][] fAValues, fABValues;
    protected final int nPtsTabulated;//tabulate fA,fAB,fC values for all graphs up to this size
    protected static final boolean checkme = false;
    protected final boolean doStatistics = false;
    public final FrequencyCounter[] sigCounter;
    long maxA, maxB, maxAB, maxC, maxQ;

    public ClusterWheatleyPartitionScreening(int nPoints, MayerFunction f) {
        this(nPoints, f, 2);
    }

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
        
        int partitionCount = 0;
        int[] iPartitions = new int[10];//start with arbitrary length and resize as needed

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
            partitionCount++;
            if(iPartitions.length < partitionCount) {
                iPartitions = Arrays.copyOf(iPartitions, partitionCount+10);
            }
            iPartitions[partitionCount-1] = j;            
        }
        if(iPartitions.length > partitionCount) {
            iPartitions = Arrays.copyOf(iPartitions, partitionCount);
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
        
    protected int[] computeAPartitions(int vs1, int i) {
        
        if((vs1 & i)==0 || (nPts[i]<3)) return null;
        
        int partitionCount = 0;
        int[] iPartitions = new int[10];//start with arbitrary length and resize as needed
        
        int iLowBit = (i&-i);//lowest bit in i
        int ii = i ^ iLowBit;//strip off low bit
        int iLow2Bit = (ii & -ii);//second-lowest bit in i
        
        if (iLowBit != vs1 && iLow2Bit != vs1) {
            //v is not in the lowest 2 bits
            
            int jBits = iLowBit | vs1;// jBits is the lowest bit and v

            //increment by the 2nd-lowest bit to ensure lowest bit is always there
            int jInc = iLow2Bit;

            //at this point jBits has (lowest bit + v)
            for (int j=jBits; j<i; j+=jInc) {//sum over partitions of i
                if ((j & jBits) != jBits) {//ensure jBits are in j
                    j |= vs1;//add v back in (v must be the missing bit, because jInc ensures lowbit is in j)
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
                partitionCount++;
                if(iPartitions.length < partitionCount) {
                    iPartitions = Arrays.copyOf(iPartitions, partitionCount+10);
                }
                iPartitions[partitionCount-1] = j;
            }
        } 
        else {
            //v is one of the lowest 2 bits
            
            // we can start at jBits and increment by the 3rd lowest bit
            int jBits = iLowBit | iLow2Bit;// jBits is the lowest 2 bits, which includes v
//            if (jBits == i) return null; // no bits left for jComp (**don't need this, because of nPts check above)
            
            int iii = ii ^ iLow2Bit;//i with 2 lowest bits off
            int jInc = (iii & -iii);//3rd lowest bit, also increment for j

            //at this point jBits has (lowest bit + v) or (v + next lowest bit)
            for (int j=jBits; j<i; j+=jInc) {//sum over partitions of i containing jBits
                // start=jBits and jInc ensure that every set includes jBits
                int jComp = i & ~j;//subset of i complementing j
                while ((j|jComp) != i && j<i) {
                    int jHighBits = j^jBits;
                    int jlow = jHighBits & -jHighBits;
                    j += jlow;
                    jComp = (i & ~j);
                }
                if (j==i) break;
                partitionCount++;
                if(iPartitions.length < partitionCount) {
                    iPartitions = Arrays.copyOf(iPartitions, partitionCount+10);
                }
                iPartitions[partitionCount-1] = j;
            }
        }
        
        if(iPartitions.length > partitionCount) {
            iPartitions = Arrays.copyOf(iPartitions, partitionCount);
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

    public ClusterAbstract makeCopy() {
        ClusterWheatleyPartitionScreening c = new ClusterWheatleyPartitionScreening(n, f);
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

    public long getCPairID() {
        return cPairID;
    }
    
    public long getLastCPairID() {
        return lastCPairID;
    }

    public double getLastValue() {
        return lastValue;
    }

    /**
     * This calculates all fQ values given that the entries for pairs have
     * already been populated.
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
     * Compute signatures used for lookup tables of fA, fAB, fC
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
                fQ[(1<<i)|(1<<j)] = (f.f(aPairs.getAPair(i,j),cPairs.getr2(i,j), beta) == -1);
            }
        }
    }
    
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
    
    //for development purposes.  Configure bonds by hand.
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

    public void setTemperature(double temperature) {
        beta = 1/temperature;
    }
    
    //class used to track the frequency of occurrence of graph signatures
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
        public void print() {
            System.out.println("Signature frequencies for n = "+n);
            double[] fractions = getFractions();
            for(int i=0; i<nfreq; i++) {
                if(freq[i]>1) System.out.println(String.format("%10d %7.5f %d", sigs[i], fractions[i], freq[i]));
            }
            System.out.println();
        }
        public double fractionTabulated() {
            int sum = 0;
            for(int i=0; i<nfreq; i++) {
                sum += freq[i];
            }
            return (double)sum/counter;
        }
        public int getn() {
            return n;
        }
        public int getEntries() {
            return entries;
        }
    }
    
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

