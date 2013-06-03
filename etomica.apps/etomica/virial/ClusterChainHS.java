package etomica.virial;



/**
 * This class calculates the sum of all chain clusters. Can be configured at construction to do rings instead.
 * 
 * @author David Kofke and Andrew Schultz 
 */
public class ClusterChainHS implements ClusterAbstract {

    protected final int n, nf;
    protected final MayerFunction f;
    protected final boolean doRing;
    
    protected final long[][] nC;
    protected final double[][] fValues;
    protected int cPairID = -1, lastCPairID = -1;
    protected long value, lastValue;
    protected double beta;
    public final boolean old = true;
    
    public ClusterChainHS(int nPoints, MayerFunction f) {
        this(nPoints, f, false);
    }
    
    public ClusterChainHS(int nPoints, MayerFunction f, boolean doRing) {
        this.n = nPoints;
        this.f = f;
        this.doRing = doRing;
        nf = 1<<n;  // 2^n

        nC = new long[n][nf];
        fValues = new double[n][n];
    }

    public ClusterAbstract makeCopy() {
        ClusterChainHS c = new ClusterChainHS(n, f);
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
        
        calcValue();
        if (Double.isNaN(value) || Double.isInfinite(value)) {
            updateF(box);
            calcValue();
            throw new RuntimeException("oops");
        }
        return value;
    }

    public long numDiagrams() {
        long savedValue = value;
        for (int i=0; i<n; i++) {
            for (int j=i+1; j<n; j++) {
                 fValues[i][j] = 1;
                 fValues[j][i] = 1;
            }
        }
        
        calcValue();
        long num = (int)Math.round(value);
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
//        if(false) {//old
//        for(int i=3; i<nf-1; i+=2) {//1-bit is always nonzero in i
//            for(int m=1; m<n; m++) {//loop over indices not in i
//                int im = 1<<m;
//                if((im & i) != 0) continue;//skip if m is in i
//                int index = i|im;
//                nC[m][index] = 0;
//                for(int k=1; k<n; k++) {//loop over indices in i
//                    int ik = 1<<k;
//                    if(ik > i) break;
////                    sum1++;
//                    if((ik & i) == 0) continue;//skip if k is not in i
////                    sum2++;
//                    nC[m][index] += fValues[m][k]*nC[k][i];
//                }
//            }
//        }
//        }
//        else {//new
                
        for(int i=3; i<nf-1; i+=2) {//1-bit is always nonzero in i
            //loop over bits not in i; start with complement of i, and in each iteration
            //get lowest bit (im=(iC&-iC)) and strip it from complement (iC^=im) until complement is empty (=0)
            for(int iC=i^(nf-1), im=(iC&-iC); iC>0; iC^=im,im=(iC&-iC)) {
                int m = log2(im);
                int iim = i|im;
                nC[m][iim] = 0;
             //loop over bits in i, in same manner as loop over complement
                for(int it=i-1, ik=(it&-it); ik>0; it^=ik,ik=(it&-it)) {
                    int k = log2(ik);
                    nC[m][iim] += fValues[m][k]*nC[k][i];
                }
            }// while(iC > 0);
        }//end for i
//        }
        
        value = 0;
        if(doRing) {
            for(int m=1; m<n; m++) {
                value += nC[m][nf-1] * fValues[m][0];
            }
            
        } else { //chains
        
            //Sum chains in which first vertex is not a leaf.
            //Consider all partitions, counting paths beginning in one partition and ending in its complement
//            if(true) {//new
//            for(int iS=3; iS<nf; iS+=4) {//keep 1 and 2 in iS-partition to prevent double counting
//                int iSComp = (nf-1)^iS;
//                for(int m=1; m<n; m++) {
//                    if(((1<<m)&iS) == 0) { //m is in iSComp
//                        for(int k=1; k<m; k++) {
//                            if(((1<<k)&iSComp) == 0) {//k is in iS
//                                value += nC[m][iSComp|1] * nC[k][iS];
//                            }
//                        }                                                
//                    } else { //m is in iS
//                        for(int k=2; k<m; k++) {
//                            if(((1<<k)&iS) == 0) { //k is in iSComp
//                                value += nC[m][iS] * nC[k][iSComp|1];
//                            }
//                        }                        
//                    }
//                }
//            }
            
            for(int iS=3; iS<nf; iS+=4) {//keep 1 and 2 in iS-partition to prevent double counting
                //loop over bits not in iS
                int iSComp = iS^(nf-1);
                for(int iC=iSComp, im=(iC&-iC); iC>0; iC^=im,im=(iC&-iC)) {
                    int m = log2(im);
                    //loop over bits in iS, in same manner as loop over complement
                    for(int it=iS-1, ik=(it&-it); ik>0; it^=ik,ik=(it&-it)) {
                        int k = log2(ik);
                        value += nC[m][iSComp|1] * nC[k][iS];
                    }
                }
            }
            
            

//            } else {//old
//            
//            for(int iS=3; iS<nf; iS+=4) {//keep 1 and 2 in i-partition to prevent double counting
//                int iSComp = (nf-1)^iS;
//                for(int m=1; m<n; m++) {
//                    if(((1<<m)&iS) == 0) continue;//skip if m is not in iS
//                    for(int k=2; k<n; k++) {
//                        if(((1<<k)&iSComp) == 0) continue;//skip if k is not in iSComp
//                        value += nC[m][iS] * nC[k][iSComp|1];
//                    }
//                }
//            }
//            }
            
            //Sum chains where first vertex is a leaf
            for(int m=1; m<n; m++) { 
                value += nC[m][nf-1];
            }
        
        }
//        System.out.println(sum1+"\t"+sum2+"\t"+(float)(sum1-sum2)/sum1+"\t"+sum3+"\t"+(float)(sum3-sum4)/sum3+"\t"+sum5+"\t"+(float)(sum5-sum6)/sum5);
        
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
        default: throw new IllegalArgumentException("Unexpected argument to log2:"+i);
        }
    }
        
    protected void updateF(BoxCluster box) {
        CoordinatePairSet cPairs = box.getCPairSet();
        AtomPairSet aPairs = box.getAPairSet();

        f.setBox(box);
        // recalculate all f values for all pairs
        for(int i=0; i<n-1; i++) {
            for(int j=i+1; j<n; j++) {
                fValues[i][j] = f.f(aPairs.getAPair(i,j),cPairs.getr2(i,j), beta);
                fValues[j][i] = fValues[i][j];
            }
        }
    }
        
    public void setTemperature(double temperature) {
        beta = 1/temperature;
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
