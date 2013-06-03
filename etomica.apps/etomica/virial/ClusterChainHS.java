package etomica.virial;



/**
 * This class calculates the sum of all chain and/or clusters for hard potential, for which the Mayer function
 * can take values of only -1 or 0. Can return a general weighted sum of ring and chain values, special-casing
 * to chain-only and ring-only values if directed at construction.
 * 
 * @author David Kofke and Andrew Schultz 
 */
public class ClusterChainHS implements ClusterAbstract {

    protected final int n, nf;
    protected final MayerFunction f;
    
    protected final long[][] nC;
    protected final double[][] fValues;
    protected final double ringFrac, chainFrac;
    protected int cPairID = -1, lastCPairID = -1;
    protected double value, lastValue;
    public final boolean old = true;
    
    /**
     * Constructs with default to perform chain-only calculation, with chainFrac = 1.0
     */
    public ClusterChainHS(int nPoints, MayerFunction f) {
        this(nPoints, f, false);
    }

    /**
     * Constructs to perform either chain-only or ring-only calculation, as directed by doRing (true for ring-only)
     */
    public ClusterChainHS(int nPoints, MayerFunction f, boolean doRing) {
        this(nPoints, f, doRing?0:1, doRing?1:0);
    }
    
    /**
     * Constructs to perform linear combination of chain and ring values, with computed cluster values
     * each weighted by given fractions, and summed to get total value
     */
    public ClusterChainHS(int nPoints, MayerFunction f, double chainFrac, double ringFrac) {
        this.n = nPoints;
        this.f = f;
        this.chainFrac = chainFrac;
        this.ringFrac = ringFrac;
        nf = 1<<n;  // 2^n

        nC = new long[n][nf];
        fValues = new double[n][n];
    }

    public ClusterAbstract makeCopy() {
        ClusterChainHS c = new ClusterChainHS(n, f);
        c.setTemperature(1);
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
        
        calcValue();
        if (Double.isNaN(value) || Double.isInfinite(value)) {
            updateF(box);
            calcValue();
            throw new RuntimeException("oops");
        }
        return value;
    }

    public long numDiagrams() {
        double savedValue = value;
        for (int i=0; i<n; i++) {
            for (int j=i+1; j<n; j++) {
                 fValues[i][j] = 1;
                 fValues[j][i] = 1;
            }
        }
        
        calcValue();
        long num = Math.round(value);
        value = savedValue;
        return num;
    }


    /*
     * Computation of sum of chain or ring diagrams.
     */
    protected void calcValue() {
        
        //nC[m][i] is the number of chains beginning at 0 and ending at m, traversing all points in i
        
        //Start with all pairwise paths from 0 to each vertex
        for(int m=1; m<n; m++) {
            nC[m][(1<<m)|1] = (int)fValues[m][0];
        }
        
        //All other paths
        //(could probably reduce memory by not including redundant first bit)
                
        for(int i=3; i<nf-1; i+=2) {//1-bit is always nonzero in i
            //the following two loops generate all pairs formed by each bit in i with each bit not in i
            
            //loop over bits not in i; start with full complement of i (i^(nf-1)), and in each iteration
            //get lowest bit (im=(iC&-iC)) and strip it from complement (iC^=im) until complement is empty (iC=0)
            for(int iC=i^(nf-1), im=(iC&-iC); iC>0; iC^=im,im=(iC&-iC)) {
                int m = log2(im);
                int iim = i|im;
                nC[m][iim] = 0;
             //loop over bits in i, in same manner as loop over complement
                for(int it=i-1, ik=(it&-it); ik>0; it^=ik,ik=(it&-it)) {
                    int k = log2(ik);
                    nC[m][iim] += fValues[m][k]*nC[k][i];
                }
            }//end for(iC)
        }//end for(i)
        
        value = 0;
        long ringValue = 0;
        long chainValue = 0;
        
        
        if(ringFrac != 0.0) {
            for(int m=1; m<n; m++) {
                ringValue += nC[m][nf-1] * fValues[m][0];
            }
        } 

        if(chainFrac != 0.0) {
        
            //Sum chains in which first vertex is not a leaf.
            //Consider all partitions, counting paths beginning in one partition and ending in its complement 
            //Use same looping structure as employed above
            for(int iS=3; iS<nf; iS+=4) {//keep 1 and 2 in iS-partition to prevent double counting
                //loop over bits not in iS
                int iSComp = iS^(nf-1);
                for(int iC=iSComp, im=(iC&-iC); iC>0; iC^=im,im=(iC&-iC)) {
                    int m = log2(im);
                    //loop over bits in iS
                    for(int it=iS-1, ik=(it&-it); ik>0; it^=ik,ik=(it&-it)) {
                        int k = log2(ik);
                        chainValue += nC[m][iSComp|1] * nC[k][iS];
                    }
                }
            }
            
            //Sum chains where first vertex is a leaf
            for(int m=1; m<n; m++) { 
                chainValue += nC[m][nf-1];
            }
            
        }//end if(chainFrac)
        
        value = chainFrac*chainValue + ringFrac*ringValue;
    }
    
    //gives position of bit for an integer having only one nonzero bit
    private final int log2(int i) {
        switch(i) {
        case (1<<0): return 0;
        case (1<<1): return 1;
        case (1<<2): return 2;
        case (1<<3): return 3;
        case (1<<4): return 4;
        case (1<<5): return 5;
        case (1<<6): return 6;
        case (1<<7): return 7;
        case (1<<8): return 8;
        case (1<<9): return 9;
        case (1<<10): return 10;
        case (1<<11): return 11;
        case (1<<12): return 12;
        case (1<<13): return 13;
        case (1<<14): return 14;
        case (1<<15): return 15;
        default: throw new IllegalArgumentException("Unexpected argument to log2: "+i);
        }
    }
        
    protected void updateF(BoxCluster box) {
        CoordinatePairSet cPairs = box.getCPairSet();
        AtomPairSet aPairs = box.getAPairSet();

        f.setBox(box);
        // recalculate all f values for all pairs
        for(int i=0; i<n-1; i++) {
            for(int j=i+1; j<n; j++) {
                fValues[i][j] = f.f(aPairs.getAPair(i,j),cPairs.getr2(i,j), 1);
                fValues[j][i] = fValues[i][j];
            }
        }
    }
        
    public void setTemperature(double temperature) {
        // we don't need no stinkin temperature!
    }

    public static void main(String[] args) {

        for(int n=5; n<13; n++) {
//            ClusterChainWheatley cc = new ClusterChainWheatley(n, null);
            ClusterSinglyConnected cs = new ClusterSinglyConnected(n, null);
            ClusterChainHS cr = new ClusterChainHS(n, null, true);
            ClusterChainHS cc2 = new ClusterChainHS(n, null);
//            cc2.old = true;
            System.out.println(n+"\t"+cs.numDiagrams()+
                    "\t"+cr.numDiagrams()+"\t"+cc2.numDiagrams());
//            cc2.old = false;
            System.out.println(n+"\t"+cs.numDiagrams()+
                    "\t"+cr.numDiagrams()+"\t"+cc2.numDiagrams());
            
        }
    }
}
