package etomica.virial;

import etomica.math.SpecialFunctions;

/**
 * This class calculates the sum of all purely singly-connected clusters using an adaptation of Wheatley's
 * recursive formulation.
 * 
 * @author David Kofke and Andrew Schultz 
 */
public class ClusterSinglyConnected implements ClusterAbstract {

    protected final int n, nf;
    protected final MayerFunction f;
    
    protected final double[][] fValues;
    protected final double[] fL, fN;
    protected final double[] bSum;
    protected int cPairID = -1, lastCPairID = -1;
    protected double value, lastValue;
    protected double beta;
    
    public ClusterSinglyConnected(int nPoints, MayerFunction f) {
        this.n = nPoints;
        this.f = f;
        fValues = new double[nPoints][nPoints];
        nf = 1<<n;  // 2^n
        fL = new double[nf];
        fN = new double[nf];
        bSum = new double[nf];
        for(int i=0; i<n; i++) {
            fL[1<<i] = 1.0;
            fN[1<<i] = 0.0;
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
     * Computation of sum of purely singly-connected diagrams.
     */
    protected void calcValue() {
        
        //Compute the fL and fN's
        // fL[i] is sum of all graphs in which high-bit node is a leaf
        // fN[i] is sum of all graphs in which high-bit node is not a leaf
        for(int m=1; m<n; m++) {//structure as nested loops so we know what high bit is
            final int iH = 1<<m; //high bit
            
            // 2 ways to compute bSum values //
            //both used stored values to increment sums as new terms are added 
            
            //compute bSum values
//            for(int i=1; i<iHighBit; i++) {
//                int j = i & -i;//lowest bit in i
//                int jj = Integer.numberOfTrailingZeros(j); // jj = log2(j)
//                if(i==j) { //only one bit
//                    bSum[m|i] = eValues[m][jj];
//                } else {
//                    int k = (i&~j);//strip off lowest bit to get previously accumulated value
//                    bSum[m|i] = bSum[m|k] + eValues[m][jj];
//                }
//            }
            
            //compute bSum values
            //in this way, loops are set up so we know what low bit is, and don't need log2
            for(int j=m-1; j>=0; j--) { //j loops down the bits to the right of m
                final int iL = 1<<j; //low bit
                final int inc = iL<<1;
                bSum[iH|iL] = fValues[m][j];
                for(int i=inc; i<iH; i+=inc) {
                    bSum[iH|i|iL] = bSum[iH|i] + bSum[iH|iL];
                }
            }
                      
            //compute fN and fL values
            for(int i=1; i<iH; i++) {
                fL[iH|i] = bSum[iH|i]*(fL[i]+fN[i]);
                fN[iH|i] = 0.0;
                final int iL = i & -i;
                final int inc = iL<<1;
                for(int iS=iL; iS<i; iS+=inc) {
                    //if ((iS & iL) == 0) continue; (enforced by loop)
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
                double x = f.f(aPairs.getAPair(i,j),cPairs.getr2(i,j), beta);
                fValues[i][j] = x;
                fValues[j][i] = x;
            }
        }
    }

    public void setTemperature(double temperature) {
        beta = 1/temperature;
    }
}
