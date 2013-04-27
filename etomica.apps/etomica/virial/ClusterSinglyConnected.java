package etomica.virial;


/**
 * This class calculates the sum of all tree clusters using an adaptation of Wheatley's
 * recursive formulation.
 * 
 * @author David Kofke and Andrew Schultz 
 */
public class ClusterSinglyConnected implements ClusterAbstract {

    protected final int n, nf;
    protected final MayerFunction f;
    
    protected final double[] fL, fN;
    protected final double[] bSum;
    protected int cPairID = -1, lastCPairID = -1;
    protected double value, lastValue;
    protected double beta;
    
    public ClusterSinglyConnected(int nPoints, MayerFunction f) {
        this.n = nPoints;
        this.f = f;
        nf = 1<<n;  // 2^n
        fL = new double[nf];
        fN = new double[nf];
        bSum = new double[nf];
        for(int i=0; i<n; i++) {
            fL[1<<i] = 1.0;
            fN[1<<i] = 0.0;
            for(int j=i+1; j<n; j++) {
                fN[(1<<i)|(1<<j)] = 0.0;
            }
        }
    }

    public ClusterAbstract makeCopy() {
        ClusterSinglyConnected c = new ClusterSinglyConnected(n, f);
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
        double savedValue = value;
//        bSum[5] = 1;
//        bSum[9] = 1;
//        bSum[10] = 1;
//        bSum[18] = 1;
//        bSum[20] = 1;
        for (int i=0; i<n; i++) {
            for (int j=i+1; j<n; j++) {
                bSum[(1<<i)|(1<<j)] = 1;
                fL[(1<<i)|(1<<j)] = bSum[(1<<i)|(1<<j)];
            }
        }
        
        calcValue();
        long num = (int)Math.round(value);
        value = savedValue;
        return num;
    }

    /*
     * Computation of sum of tree diagrams.
     */
    protected void calcValue() {
        
        //Compute the fL and fN's
        // fL[i] is sum of all graphs in which low-bit node is a leaf
        // fN[i] is sum of all graphs in which low-bit node is not a leaf
        for(int i=0; i<nf; i++) fN[i] = 0;
        for(int m=2; m<n; m++) {//structure as nested loops so we know what high bit is
            final int iH = 1<<m; //high bit
            
            //bSum[i] is the sum of all bonds formed from the low-bit of i with each other non-zero bit of i
            
            //calculation of bSum is performed for any i by adding the high-low (iH|iL) bit interaction to the sum
            //obtained without iL, obtained from a previous iteration
            for(int i=iH+3; i<(iH<<1); i++) {
                int iL = i & -i;//low bit
                int i0 = i^iL;//i, without the low bit
                int iH0 = i^iH;
                if(i0 == iH) continue;//only two bits in i; we skip this because we start with all pairs in bSum and fL
                bSum[i] = bSum[iH0] + bSum[iH|iL];

                //compute fN and fL values
                fL[i] = bSum[i]*(fL[i0]+fN[i0]);
                //fN[i] = 0.0;
                int inc = i0 & -i0;
                for(int iS=iH+iL; iS<i; iS+=inc) {//structure loop to force iS to contain iH and iL bits
                    int iSComp = i & ~iS;
                    if ((iSComp | iS) != i) continue;
                    fN[i] += fL[iS]*(fL[iL|iSComp] + fN[iL|iSComp]);
                }
            }
        }

        value = fL[nf-1] + fN[nf-1];
    }

    protected void updateF(BoxCluster box) {
        CoordinatePairSet cPairs = box.getCPairSet();
        AtomPairSet aPairs = box.getAPairSet();

        f.setBox(box);
        // recalculate all f values for all pairs
        for(int i=0; i<n-1; i++) {
            for(int j=i+1; j<n; j++) {
                int index = (1<<i)|(1<<j);
                bSum[index]= f.f(aPairs.getAPair(i,j),cPairs.getr2(i,j), beta);
                fL[index] = bSum[index];
            }
        }
    }

    public void setTemperature(double temperature) {
        beta = 1/temperature;
    }
    
    public static void main(String[] args) {
        for(int i=0; i<100; i++) {
            System.out.println(Integer.toBinaryString((1<<10)|i));
        }
        ClusterChain cc = new ClusterChain(5, null);
        ClusterSinglyConnected cs = new ClusterSinglyConnected(5, null);
        System.out.println(cc.numDiagrams());
        System.out.println(cs.numDiagrams());

    }
}
