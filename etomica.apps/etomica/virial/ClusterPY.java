package etomica.virial;

import etomica.math.SpecialFunctions;

public class ClusterPY extends ClusterWheatley {

    protected final double[][][] t, h, c;

    public ClusterPY(int nPoints, MayerFunction f) {
        super(nPoints, f);
        
        int nf = 1<<n;
        h = new double[n][n][nf];
        c = new double[n][n][nf];
        t = new double[n][n][nf];
        
    }

    //compute PY diagrams and subtract from doubly-connected diagrams computed in superclass
    protected void calcValue() {
        super.calcValue();
        
        //fill in c0 and h0 with pair f-bonds
        for(int i=0; i<n-1; i++) {
            for(int j=i+1; j<n; j++) {
                int index = (1<<i)|(1<<j);
                h[j][i][index] = fQ[index]-1;
                h[i][j][index] = fQ[index]-1;
                c[j][i][index] = fQ[index]-1;
                c[i][j][index] = fQ[index]-1;
            }
        }
        
        int nf = 1<<n;
iLoop:  for(int i=1; i<nf; i++) {//sum over subsets of points

            //compute tn = sum[ c_j h_{n-j-1}, j=0,n-1]
            for(int iS=1; iS<i; iS++) {//sum over partitions of i
                int iSComp = i & ~iS;
                if ((iSComp | iS) != i) continue;

                for(int jL=0; jL<n; jL++) {//leaf on one partition
                    int iL = 1<<jL;
                    if((iL & iS) == 0) continue;
                    for(int jR=0; jR<n; jR++) {//leaf on the other partition
                        int iR = 1<<jR;
                        if ((iR|iL) == i) continue iLoop;
                        if((iR & iSComp) == 0) continue;
                        t[jL][jR][i] = 0.0;
                        for(int jM=0; jM<n; jM++) {//leaf where chains are spliced
                            int iM = 1<<jM;
                            if(jM==jL || jM==jR || (iM&i)==0) continue;
                            t[jL][jR][i] += c[jL][jM][iS|iM]*h[jM][jR][iSComp|iM];
                        }
                    }
                }
            }
            //PY: c_n = f * t_n
            for(int jL=0; jL<n; jL++) {
                for(int jR=0; jR<n; jR++) {
                    t[jR][jL][i] = t[jL][jR][i];
                    c[jL][jR][i] = t[jL][jR][i] * (fQ[(1<<jR)|(1<<jL)] - 1); 
                    c[jR][jL][i] = c[jL][jR][i];
                    h[jL][jR][i] = t[jL][jR][i] + c[jL][jR][i]; 
                    h[jR][jL][i] = h[jL][jR][i];
                }
            }
        }
        
        double sum = 0.0;
        for(int jR=1; jR<n; jR++) {
            for(int jL=0; jL<jR; jL++) {
                sum += c[jR][jL][nf-1];
            }
        }

        value -= -sum/n;

    }

    
}
