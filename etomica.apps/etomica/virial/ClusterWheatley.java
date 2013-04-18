package etomica.virial;

import java.util.HashSet;
import java.util.List;

import etomica.graph.model.Bitmap;
import etomica.graph.model.BitmapFactory;
import etomica.graph.model.Coefficient;
import etomica.graph.model.Edge;
import etomica.graph.model.EdgeVisitor;
import etomica.graph.model.Graph;
import etomica.graph.model.Node;
import etomica.graph.model.NodeVisitor;
import etomica.graph.property.IsBiconnected;

/**
 * This class calculates the sum of all biconnected clusters using Wheatley's
 * recursive formulation.
 * 
 * @author David Kofke and Andrew Schultz 
 */
public class ClusterWheatley implements ClusterAbstract {

    protected final int n;
    protected final MayerFunction f;
    
    protected final double[] fQ, fC;
    protected final double[] fA, fB;
    protected int cPairID = -1, lastCPairID = -1;
    protected double value, lastValue;
    protected double beta;
    protected final Bitmap bondMap;
    protected final MyGraph myGraph;
    protected final IsBiconnected isBi;
    protected boolean doBiconCheck, doCliqueCheck;
    protected final byte[] outDegree;

    public ClusterWheatley(int nPoints, MayerFunction f) {
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
        bondMap = BitmapFactory.createBitmap((byte)n, false);
        outDegree = new byte[n];
        myGraph = new MyGraph((byte)nPoints, bondMap, outDegree);
        isBi = new IsBiconnected();
        doBiconCheck = n>2;
        doCliqueCheck = n>3;
    }
    
    public void setDoBiconCheck(boolean newDoBiconCheck) {
        doBiconCheck = newDoBiconCheck && n>2;
    }

    public void setDoCliqueCheck(boolean newDoCliqueCheck) {
        doCliqueCheck = newDoCliqueCheck && n>3;
    }

    public ClusterAbstract makeCopy() {
        ClusterWheatley c = new ClusterWheatley(n, f);
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
            //loop over pairs formed from j and each point in partition; multiply by bond for each pair
            //all such pairs will be with bits higher than j, as j is the lowest bit in i
            for (int l=(j<<1); l<i; l=(l<<1)) {
                if ((l&i)==0) continue; //l is not in partition
                fQ[i] *= fQ[l | j];
            }
        }
    }

    /*
     * Computation of sum of biconnected diagrams.
     */
    protected void calcValue(BoxCluster box) {
        int nf = 1<<n;
        if (false && total>0 && total%1000000==0) {
            System.out.println("screened: "+((double)screened)/total+"     zero-not-screened: "+((double)(total-screened-notzero))/(total-screened));
        }
        total++;
        screened++;
        edgeCount = 0;
        if (doBiconCheck) {
            int bit=0;
            for (int i=0; i<n; i++) {
                outDegree[i] = 0;
            }

            for (int i=0; i<n-1; i++) {
                for (int j=i+1; j<n; j++) {
                    boolean fBond = (fQ[(1<<i)|(1<<j)] == 0); 
                    if (fBond) {
                        bondMap.setBit(bit);
                        outDegree[i]++;
                        outDegree[j]++;
                        edgeCount++;
                    }
                    else {
                        bondMap.clearBit(bit);
                    }
                    bit++;
                }
            }
            for (int i=0; i<n; i++) {
                if (outDegree[i] < 2) {
                    value = 0;
                    return;
                }
            }

            if (doCliqueCheck) {
                for (byte i=0; i<n; i++) {
                    if (outDegree[i] < n-1) {
                        boolean isClique = true;
jLoop:                  for (byte j=0; j<outDegree[i]-1; j++) {
                            byte jj = myGraph.getOutNode(i, j);
                            if (outDegree[jj] < outDegree[i]-1) continue;
                            for (byte k=(byte)(j+1); k<outDegree[i]; k++) {
                                byte kk = myGraph.getOutNode(i, k);
                                if ((fQ[(1<<kk)|(1<<jj)] != 0)) {
                                    isClique = false;
                                    break jLoop;
                                }
                            }
                        }
                        if (isClique) {
                            value = 0;
                            return;
                        }
                    }
                }
            }
        }
        screened--;
        calcFullFQ(box);

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
                        if (j==i) break;
                        fA[i] += fB[j] * (fB[jComp|vs1] + fA[jComp|vs1]);
                    }
                }

                fB[i] -= fA[i];//remove from B graphs that contain articulation point at v
            }
        }

        value = (1-n)*fB[nf-1]; ///SpecialFunctions.factorial(n);
        if (value != 0) {
            notzero++;
        }
    }
    
    /**
     * Returns edgeCount (number of overlaps) of configuration passed to
     * checkConfig
     */
    public int getEdgeCount() {
        return edgeCount;
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

    /**
     * This is a minimal implemenation of Graph that is sufficient to handle
     * the biconnectivity test (via IsBiconnected)
     * 
     * @author Andrew Schultz
     */
    public static class MyGraph implements Graph {

        protected final Bitmap bitmapStore;
        protected final byte n;
        protected final byte[] outDegree;

        /*
         * We just wrap the bitmapStore and outDegree.
         * ClusterWheatley will update both as needed.
         */
        public MyGraph(byte n, Bitmap bitmapStore, byte[] outDegree) {
            this.n = n;
            this.bitmapStore = bitmapStore;
            this.outDegree = outDegree;
        }

        public int compareTo(Graph o) {throw new RuntimeException("don't be calling me");}
        public Coefficient coefficient() {throw new RuntimeException("don't be calling me");}
        public int[] factors() {throw new RuntimeException("don't be calling me");}
        public void addFactors(int[] newFactors) {throw new RuntimeException("don't be calling me");}
        public void setNumFactors(int numFactors) {throw new RuntimeException("don't be calling me");}
        public Graph copy() {throw new RuntimeException("don't be calling me");}
        public void deleteEdge(byte edgeId) {throw new RuntimeException("don't be calling me");}
        public void deleteEdge(byte fromNode, byte toNode) {throw new RuntimeException("don't be calling me");}
        public List<Edge> edges() {throw new RuntimeException("don't be calling me");}
        public String edgesToString() {throw new RuntimeException("don't be calling me");}
        public Edge getEdge(byte edgeId) {throw new RuntimeException("don't be calling me");}
        public Edge getEdge(byte fromNode, byte toNode) {throw new RuntimeException("don't be calling me");}

        public byte edgeCount() {
            return (byte)bitmapStore.bitCount();
        }

        public byte getEdgeId(byte fromNode, byte toNode) {
            byte id0 = fromNode;
            byte id1 = toNode;
            if (id0 > id1) {
                id0 = toNode;
                id1 = fromNode;
            }
            return (byte)((2*n-id0-1)*id0/2 + (id1-id0-1));
        }

        public byte getFromNode(byte edgeId) {
            for (int i=0; i<n; i++) {
                if ((2*n-i-1)*i/2 > edgeId) {
                    return (byte)(i-1);
                }
            }
            throw new RuntimeException("invalid edgeID");
        }

        public byte getOutDegree(byte node) {
            // called
            return outDegree[node];
        }

        public byte getOutNode(byte node, byte index) {
            // called
            byte c = 0;
            int bit = node-1;
            for (byte i=0; i<node; i++) {
                if (bitmapStore.testBit(bit)) {
                    if (c==index) return i;
                    c++;
                }
                bit += n-i-2;
            }
            bit++;
            for (byte i=(byte)(node+1); i<n; i++) {
                if (bitmapStore.testBit(bit)) {
                    if (c==index) return i;
                    c++;
                }
                bit++;
            }
            throw new RuntimeException("not that many bonds");
        }

        public Bitmap getStore() {
            return bitmapStore;
        }

        public byte getToNode(byte edgeId) {
            for (int i=0; i<n; i++) {
                int foo = (2*n-i-1)*i/2;
                if (foo > edgeId) {
                    return (byte)(i+(foo-edgeId));
                }
            }
            throw new RuntimeException("invalid edgeID");
        }

        public boolean hasEdge(byte edgeId) {
            return bitmapStore.testBit(edgeId);
        }

        public boolean hasEdge(byte fromNode, byte toNode) {
            return hasEdge(getEdgeId(fromNode, toNode));
        }

        public byte nodeCount() {
            // called
            return n;
        }

        public Node getNode(byte node) {throw new RuntimeException("don't be calling me");}
        public String getSignature() {throw new RuntimeException("don't be calling me");}
        public List<Node> nodes() {throw new RuntimeException("don't be calling me");}
        public String nodesToString() {throw new RuntimeException("don't be calling me");}
        public void putEdge(byte edgeId) {throw new RuntimeException("don't be calling me");}
        public void putEdge(byte fromNode, byte toNode) {throw new RuntimeException("don't be calling me");}
        public String toSVG(int dim) {throw new RuntimeException("don't be calling me");}
        public void visitEdges(EdgeVisitor visitor) {throw new RuntimeException("don't be calling me");}
        public void visitNodes(NodeVisitor visitor) {throw new RuntimeException("don't be calling me");}
        public void createReverseEdges() {throw new RuntimeException("don't be calling me");}
    }
}
