package etomica.virial;

import java.util.BitSet;
import java.util.HashSet;

import etomica.graph.model.Graph;
import etomica.graph.model.impl.GraphImpl;
import etomica.graph.operations.MaxIsomorph;
import etomica.graph.operations.GraphOp.GraphOpNull;
import etomica.graph.operations.MaxIsomorph.MaxIsomorphParameters;

/**
 * This class calculates the sum of all biconnected clusters using Wheatley's
 * recursive formulation.
 * 
 * @author David Kofke and Andrew Schultz 
 */
public class ClusterWheatleyPartitionScreening implements ClusterAbstract {

    protected final int n, nf;
    protected final MayerFunction f;
    
    protected final double[] fQ, fC;
    protected final double[] fA, fB;
    protected int cPairID = -1, lastCPairID = -1;
    protected double value, lastValue;
    protected double beta;
    protected final byte[] outDegree;
    protected final int[] fullBondMask;
    protected final boolean[] cliqueSet;
    protected int cliqueCount, eCliqueCount;
    protected boolean precalcQ = true;
    protected final int[] cliqueList;
    
    protected final int[] nBonds;
    protected final int[] nPts;
    protected final int[] vCount;
    protected final int[] sig;
    long bigSum1, bigSum2, count, bigSumN, bigSumNCount;
    public final static int[] fAVals = new int[13];
    public final static int[] fBVals = new int[13];

    protected final int[] zero = {0,0,0,0,0,0,0,0,0,0,0,0,0};
    protected final int[][][] fAList;
//                                   = {{0},                      //0
//                                      {0},                      //1
//                                      {0},                      //2
//                                      {0,+1, 1, 1},             //3
//                                      {0,-4,-2, 0, 2},          //4
//                                      {0,+18, 6, 0, 0, 6},      //5
//                                      {0,-96,-24,0,0,0,24},     //6
//                                      {0,+600,120,0,0,0,0,120},  //7
//                                      {0,-4320,-720,0,0,0,0,0,720},//8
//                                      {0,+35280,5040,0,0,0,0,0,0,5040},//9
//                                      {0, -322560, -40320, 0, 0, 0, 0, 0, 0, 0, 40320},//10
//                                      {0, 3265920, 362880, 0, 0, 0, 0, 0, 0, 0, 0, 362880},//11
//                                      {0, -36288000, -3628800, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3628800},//12
//                                      {0,+439084800,39916800,0,0,0,0,0,0,0,0,0,0,39916800}//13
//                                      };
    protected final int[][][] fBList; 
//                                   = {{0},                      //0
//                                      {0},                      //1
//                                      {0},                      //2
//                                      {+2, 1, 0,-1},            //3
//                                      {-6,-2, 0, 0,-2},         //4
//                                      {+24, 6, 0, 0, 0,-6},      //5
//                                      {-120,-24,0,0,0,0,-24},    //6
//                                      {+720,+120,0,0,0,0,0,-120},//7
//                                      {-5040,-720,0,0,0,0,0,0,-720},//8
//                                      {+40320,+5040,0,0,0,0,0,0,0,-5040},//9
//                                      {-362880, -40320, 0, 0, 0, 0, 0, 0, 0, 0, -40320},//10
//                                      {3628800, 362880, 0, 0, 0, 0, 0, 0, 0, 0, 0, -362880},//11
//                                      {-39916800, -3628800, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -3628800},//12    
//                                      {479001600,+39916800,0,0,0,0,0,0,0,0,0,0,-39916800}//13
//                                      };
    protected final int[][] fCList;
    
    protected final int[][] fAValues, fBValues;
    protected final int nPtsTabulated = 5;//tabulate fA,fB,fC values for all graphs up to this size
    protected static final boolean checkme = false;
    protected final boolean doStatistics = false;
    public final FrequencyCounter[] sigCounter;

    public ClusterWheatleyPartitionScreening(int nPoints, MayerFunction f) {
        this.n = nPoints;
        this.f = f;
        nf = 1<<n;  // 2^n
        fQ = new double[nf];
        fC = new double[nf];
        nBonds = new int[nf];
        nPts = new int[nf];
        vCount = new int[nf];
        fAValues = new int[nf][];
        fBValues = new int[nf][];
        sig = new int[nf];

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
            fQ[1<<i] = 1.0;
            nBonds[1<<i] = 0;
        }
        fA = new double[nf];
        fB = new double[nf];
        outDegree = new byte[n];
        fullBondMask = new int[nf];
        cliqueSet = new boolean[nf];
        cliqueList = new int[nf];
        
        fAList = new int[n+1][][];
        fBList = new int[n+1][][];
        fCList = new int[n+1][];
        for(int np=2; np<=nPtsTabulated; np++) {
            fAList[np] = new int[1<<nPairs(np)][np];
            fBList[np] = new int[1<<nPairs(np)][np];
            fCList[np] = new int[1<<nPairs(np)];
        }
        int sigMax = (1 << nPairs(nPtsTabulated));
        int bond = 0;
        int noBond = 1;
        //construct table of fA, fB, and fC for each signature for each nPts being tabulated
        for(int s=1; s<sigMax; s++) {//loop over signatures
            for(int np=2; np<=nPtsTabulated; np++) {
                if(s >= (1<<nPairs(np))) continue;//signature has more bonds than can be formed by np points
                //construct bonds from signature
                int bondBit = 1;
                for(int j=(1<<(np-2)); j>0; j=(j>>1)) {
                    for (int l=(j<<1); l<(1<<np); l=(l<<1)) {
                        fQ[l | j] = ((s & bondBit) == bondBit) ? bond : noBond;
                        bondBit = bondBit << 1;
                    }
                }

                calcFullFQ(null);

//                System.out.println(np+"\t"+s+"\t"+sig[(1<<np)-1]+"\t"+(s-sig[(1<<np)-1]));
            
                calcfC(false);
                int i = (1<<np) - 1;
                fCList[np][s] = (int)fC[i];
                for(int v=0; v<np; v++) {
                    calcfAB(v, false);
                    fAList[np][s][v] = (int)fA[i];
                    fBList[np][s][v] = (int)fB[i];
                }//end of v-loop
            }//end of np-loop
        }//end of s-loop
        if(checkme) System.out.println("Note -- checkme is true");
        System.out.println("tabulation complete");
//        System.exit(1);
        
        if(doStatistics) {
            System.out.println("Note -- doStatistics is true");
            for(int i=0; i<n; i++) {
                sigCounter[i].reset();
            }
        }        
   }
    
    private int nPairs(int nPts) {
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
      int thisCPairID = cPairs.getID();
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
        eCliqueCount = 0;
        // generate all partitions and compute product of e-bonds for all pairs in partition
        
        for (int i=3; i<nf; i++) {
            vCount[i] = 0;
            if(nPts[i] == 1) continue;//1-point sets filled on construction
            if(nPts[i] == 2) {// 2-point set fQ's were filled when bonds were computed, so skip
                if(fQ[i] == 0) {
                    nBonds[i] = 1;//fQ is e-bond, so if zero it means f-bond is nonzero
                    sig[i] = 1;
                } else {
                    nBonds[i] = 0;
                    sig[i] = 0;
                }
                if(nPtsTabulated >= 2) {
                    fAValues[i] = fAList[2][sig[i]];
                    fBValues[i] = fBList[2][sig[i]];
                    fC[i] = fCList[2][sig[i]];
                }
                if(doStatistics) sigCounter[1].add(sig[i]);
                continue;               
            }
            int j = i & -i;//lowest bit in i
            int k = i&~j; //strip j bit from i and set result to k
            fQ[i] = fQ[k]; //initialize with previously-computed product of all pairs in partition, other than j
            nBonds[i] = nBonds[k];
            sig[i] = sig[k];
            //loop over pairs formed from j and each point in partition; multiply by bond for each pair
            //all such pairs will be with bits higher than j, as j is the lowest bit in i
            
            int bondBit = 0;
            boolean doSig = (nPts[i] <= nPtsTabulated); 
            if(doSig || doStatistics) bondBit = 1 << nPairs(nPts[k]);
            for (int l=(j<<1); l<i; l=(l<<1)) {
                if ((l&i)==0) continue; //l is not in partition
                fQ[i] *= fQ[l | j];
                if(fQ[l | j] == 0) {
                    nBonds[i]++;
                    if(doSig || doStatistics) sig[i] |= bondBit;
                } 
                if(doSig || doStatistics) bondBit = bondBit << 1;
            }
            if(doSig) {
                fAValues[i] = fAList[nPts[i]][sig[i]];
                fBValues[i] = fBList[nPts[i]][sig[i]];
                fC[i] = fCList[nPts[i]][sig[i]];
            }
            
            if (fQ[i] != 0) eCliqueCount++;
            
            if(doStatistics) sigCounter[nPts[i]-1].add(sig[i]);
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
        
        int sum1 = 0;
        int sum2 = 0;
        int counter = 0;

        calcfC(true);
        
        for(int i=1; i<nf; i++) {
            vCount[i] = 0;
        }

        for (int v=0; v<n; v++) {
            calcfAB(v, true);
            
           if (checkme) {
                for(int i=0; i<nf; i++) {
                    if(fB[i]==0 && nPts[i] > 2) sum2++;
                    counter++;
                }
                for(int i=0; i<nf; i++) {
                }
            }//end of checkme
        }//end of v-loop
        
        value = (1-n)*fB[nf-1];

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
            if(fBValues[i] != null & useTable) {
                continue;//fC was set in calcFullFQ
            }
            
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
    }
    
    protected final void calcfAB(int v, final boolean useTable) {
        
        if(v == 0) {
            if(checkme && useTable) {
                calcfAB0(false);
                for(int i=1; i<nf; i++) {
                    if(fAValues[i] != null) {
                        int fATable = (i%2==1) ? fAValues[i][vCount[i]] : 0;
                        int fBTable = (i%2==1) ? fBValues[i][vCount[i]++] : (int)fC[i];
                        int fAerr = (int)(fA[i]-fATable);
                        int fBerr = (int)(fB[i]-fBTable);
                        if(fAerr != 0 || fBerr != 0) System.out.println(v+"\t"+i+"\t"+(int)fA[i]+"\t"+fATable+"\t"+fAerr+"\t"+(int)fB[i]+"\t"+fBTable+"\t"+fBerr);
                    }
                }
                return;
            }
            
            calcfAB0(useTable);
            return;
        }

        int vs1 = 1<<v;
        for (int i=vs1+1; i<nf; i++) {
            
//          fB[v][i] = fB[v-1][i];//no a.p. at v or below, starts with those having no a.p. at v-1 or below
            //rest of this is to generate A (diagrams having a.p. at v but not below), and subtract it from B
            if ((i & vs1) == 0) {
                fA[i] = 0;
                continue;//if i doesn't contain v, fA and fB are done
            }
            
            if(fAValues[i] != null & useTable) {
                if(checkme) {
                    calcfAB(vs1, i);
                    int fAerr = (int)(fA[i]-fAValues[i][vCount[i]]);
                    int fBerr = (int)(fB[i]-fBValues[i][vCount[i]]);
                    if(fAerr != 0 || fBerr != 0) System.out.println(v+"\t"+i+"\t"+(int)fA[i]+"\t"+fAValues[i][vCount[i]]+"\t"+fAerr+"\t"+(int)fB[i]+"\t"+fBValues[i][vCount[i]]+"\t"+fBerr);
                }
                fA[i] = fAValues[i][vCount[i]];
                fB[i] = fBValues[i][vCount[i]];
                vCount[i]++;
            } else {
                calcfAB(vs1, i);
            }

        }//end of i-loop
    }

    //used only by calcfAB (note that this takes vs1 and other calcfAB takes v)
    protected final void calcfAB(int vs1, int i) {
        
        fA[i] = 0;
        int iLowBit = (i&-i);//lowest bit in i
        if (iLowBit == i) { //lowest bit is only bit; fA and fB are done
            return;
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
            if (jBits == i) return; // no bits left for jComp
            
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
                fA[i] += fB[j] * (fB[jComp|vs1] + fA[jComp|vs1]);
            }
        }
        fB[i] -= fA[i];//remove from B graphs that contain articulation point at v
    }
    
    protected final void calcfAB0(final boolean useTable) {

        for (int i=2; i<nf; i+=2) {
            // all even sets don't contain 1
            fA[i] = 0;
            fB[i] = fC[i];
        }

        fA[1] = 0;
        fB[1] = fC[1];
        for (int i=3; i<nf; i+=2) {
            // every set will contain 1
            if(fAValues[i] != null && useTable) {
                fA[i] = fAValues[i][0];
                fB[i] = fBValues[i][0];
                vCount[i] = 1;
            } else {
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
    
    //for development purposes.  Configure bonds by hand.
    protected void updateF() {
        for(int i=0; i<n-1; i++) {
            for(int j=i+1; j<n; j++) {
                fQ[(1<<i)|(1<<j)] = 0.0;//f.f(aPairs.getAPair(i,j),cPairs.getr2(i,j), beta)+1;
            }
        }
        fQ[(1<<2)|(1<<0)] = 1.0;
        fQ[(1<<0)|(1<<1)] = 1.0;
        fQ[(1<<1)|(1<<2)] = 1.0;
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
        cw.calcFullFQ(null);
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

