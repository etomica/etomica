package etomica.virial;

/**
 * This class calculates the sum of all chain clusters using an adaptation of Wheatley's
 * recursive formulation.
 * 
 * @author David Kofke and Andrew Schultz 
 */
public class ClusterChain implements ClusterAbstract {

    protected final int n, nf;
    protected final MayerFunction f;
    
    protected final double[][][] fL;
    protected final double[] bond;
    protected int cPairID = -1, lastCPairID = -1;
    protected double value, lastValue;
    protected double beta;
    
    public ClusterChain(int nPoints, MayerFunction f) {
        this.n = nPoints;
        this.f = f;
        nf = 1<<n;  // 2^n
        fL = new double[n][n][nf];
        bond = new double[nf];
        for(int i=0; i<n; i++) {
                fL[i][i][1<<i] = 1.0;
        }
    }

    public ClusterAbstract makeCopy() {
        ClusterChain c = new ClusterChain(n, f);
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

    /*
     * Computation of sum of purely singly-connected diagrams.
     */
    protected void calcValue() {
        
        //Compute the fL and fN's
        // fL[j][k][i] is sum of all permutations of partition i in which j and k are at the ends
        for(int m=1; m<n; m++) {//structure as nested loops so we know what high bit is
            final int iH = 1<<m; //high bit
                                              
            for(int i=1; i<iH; i++) {
                for(int j=1; j<m; j++) {
                    if()
                    fL[m][j][iH|i] = bond[iH|i]*(fL[i]+fN[i]);
                }
                final int iLowBit = i & -i;
                for(int iS=1; iS<i; iS++) {
                    if ((iS & iLowBit) == 0) continue;
                    final int iSComp = i & ~iS;
                    if ((iSComp | iS) != i) continue;
                    fN[iH|i] += fL[iH|iS]*(fL[iH|iSComp] + fN[iH|iSComp]);
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
                bond[(1<<i)|(1<<j)] = f.f(aPairs.getAPair(i,j),cPairs.getr2(i,j), beta)+1;
            }
        }
    }

    public void setTemperature(double temperature) {
        beta = 1/temperature;
    }
}
