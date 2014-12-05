/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial;

import etomica.math.SpecialFunctions;

public class ClusterPY extends ClusterWheatleyHS {

    protected final double[][][] t, h, c;

    public ClusterPY(int nPoints, MayerFunction f) {
        super(nPoints, f);
        
        int nf = 1<<n;
        h = new double[n][n][nf];
        c = new double[n][n][nf];
        t = new double[n][n][nf];
        
    }

    //compute PY diagrams and subtract from doubly-connected diagrams computed in superclass
    public double calcValue(BoxCluster box) {
        super.calcValue(box);
        value /= SpecialFunctions.factorial(n);
        
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
        for(int i=0; i<n; i++) {
            for(int j=0; j<n; j++) {
                for (int k=0; k<nf; k++) {
                    t[i][j][k] = 0;
                }
            }
        }
        
        for(int i=1; i<nf; i++) {//sum over subsets of points
            int iSize = Integer.bitCount(i);
            if (iSize<3) continue;
            //compute tn = sum[ c_j h_{n-j-1}, j=0,n-1]
            for(int iS=1; iS<i; iS++) {//sum over partitions of i
                int iSComp = i & ~iS;
                if ((iSComp | iS) != i) continue;
                int iSSize = Integer.bitCount(iS);
                if (iSSize<2) continue;

                // fac is the number of ways to choose cj hk (in the loops below)
                long fac = SpecialFunctions.factorial(iSize-2)/(SpecialFunctions.factorial(iSSize-2)*SpecialFunctions.factorial(iSize-iSSize-1));

                for(int jL=0; jL<n; jL++) {//leaf on one partition
                    int iL = 1<<jL;
                    if((iL & iS) == 0) continue;
                    for(int jR=0; jR<n; jR++) {//leaf on the other partition
                        int iR = 1<<jR;
                        if((iR & iSComp) == 0) continue;
                        for(int jM=0; jM<n; jM++) {//leaf where chains are spliced
                            int iM = 1<<jM;
                            if(jM==jL || jM==jR || (iM&iS)==0) continue;
                            t[jL][jR][i] += c[jL][jM][iS]*h[jM][jR][iSComp|iM]/fac;
                        }
                    }
                }
            }
            //PY: c_n = f * t_n
            for(int jL=0; jL<n; jL++) {
                for(int jR=0; jR<n; jR++) {
                    c[jL][jR][i] = t[jL][jR][i] * (fQ[(1<<jR)|(1<<jL)] - 1); 
                    h[jL][jR][i] = t[jL][jR][i] + c[jL][jR][i]; 
                }
            }
        }

        double sum = 0.0;
        for(int jR=0; jR<n; jR++) {
            for(int jL=0; jL<n; jL++) {
                if (jR==jL) continue;
                sum += c[jR][jL][nf-1];
            }
        }

        // we need to divide by not only (-1/n) but also the number of pairs above
        // (each pair would give us the appropriate value, we use all pairs to get
        // all permutations)
        value -= -sum/(n*n*(n-1));
        return value;
    }

    
}
