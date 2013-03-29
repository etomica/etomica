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

    public long numDiagrams() {
        for (int i=0; i<n-1; i++) {
            for (int j=i+1; j<n; j++) {
                int k = (1<<i) | (1<<j);
                fL[i][j][k] = 1;
                fL[j][i][k] = 1;
            }
        }
        double oldValue = value;
        calcValue();
        long num = (long)value;
        value = oldValue;
        return num;
    }

    /*
     * Computation of sum of purely singly-connected diagrams.
     */
    protected void calcValue() {
        
        //Compute the fL and fN's
        // fL[j][k][i] is sum of all permutations of partition i in which j and k are at the ends
iLoop:  for(int i=1; i<nf; i++) {//sum over subsets of points
            for(int iS=1; iS<i; iS++) {//sum over partitions of i
                int iSComp = i & ~iS;
                if ((iSComp | iS) != i) continue;

                for(int jL=0; jL<n; jL++) {//leaf on one partition
                    int iL = 1<<jL;
                    if((iL & iS) == 0) continue;
                    for(int jR=0; jR<n; jR++) {//leaf on the other partition
                        int iR = 1<<jR;
                        if((iR & iSComp) == 0) continue;
                        if ((iL|iR) == i) continue iLoop;
                        fL[jL][jR][i] = 0.0;
                        for(int jM=0; jM<n; jM++) {//leaf where chains are spliced
                            int iM = 1<<jM;
                            if(jM==jL || jM==jR || (iM&i)==0) continue;
                            fL[jL][jR][i] += fL[jL][jM][iS|iM]*fL[jM][jR][iSComp|iM];
                        }
                        fL[jR][jL][i] = fL[jL][jR][i];
                    }
                }
            }
        }
//        for(int m=1; m<n; m++) {//structure as nested loops so we know what high bit is
//            final int iH = 1<<m; //high bit
//                                              
//            for(int i=1; i<iH; i++) {//loop over partitions of j<m
//                final int iLowBit = i & -i;
//                final int inc = iLowBit<<1;
//                for(int iS=iLowBit; iS<i; iS+=inc) {
//                    //if ((iS & iLowBit) == 0) continue;
//                    final int iSComp = i & ~iS;
//                    if ((iSComp | iS) != i) continue;
//
//                    for(int jL=1; jL<m; jL++) {
//                        int iL = 1<<jL;
//                        if((iL & iS) == 0) continue;
//                        for(int jR=0; jR<jL; jR++) {
//                            int iR = 1<<jR;
//                            if((iR & iSComp) == 0) continue;
//                            for(int jM=1; jM<m; jM++) {
//                                int iM = 1<<jM;
//                                if(jM==jL || jM==jR || (iM&i)==0) continue;
//                                fL[jL][jR][iH|i] = fL[jL][jM][iH|]*fL[jL][jR][i];
//                                fL[jL][jR][iH|i] = fL[jL][]
//                            }
//                        }
//                    }
//                    final int iLowBit = i & -i;
//                for(int iS=1; iS<i; iS++) {
//                    if ((iS & iLowBit) == 0) continue;
//                    final int iSComp = i & ~iS;
//                    if ((iSComp | iS) != i) continue;
//                    fN[iH|i] += fL[iH|iS]*(fL[iH|iSComp] + fN[iH|iSComp]);
//                }
//            }
//        }

        //sum over all leaf pairs for full-length chain
        value = 0.0;
        for(int jR=1; jR<n; jR++) {
            for(int jL=0; jL<jR; jL++) {
                value += fL[jR][jL][nf-1];
            }
        }
    }

    protected void updateF(BoxCluster box) {
        CoordinatePairSet cPairs = box.getCPairSet();
        AtomPairSet aPairs = box.getAPairSet();

        f.setBox(box);
        // recalculate all f values for all pairs
        for(int i=0; i<n-1; i++) {
            for(int j=i+1; j<n; j++) {
                int index = (1<<i)|(1<<j);
                fL[j][i][index] = f.f(aPairs.getAPair(i,j),cPairs.getr2(i,j), beta);
                fL[i][j][index] = fL[j][i][index];
            }
        }
    }

    public void setTemperature(double temperature) {
        beta = 1/temperature;
    }
}
