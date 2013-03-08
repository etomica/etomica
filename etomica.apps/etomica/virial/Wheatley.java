package etomica.virial;

import java.util.Arrays;

/* 
 * Implementation of Wheatley's algorithm for computing sum of biconnected diagrams. For a graph of order n
 * the partitions are represented by bit strings of length n, formed as integers from 1 to 2^n. Thus, for example,
 * for n = 4, the bits 1101 (integer 13) indicates the partition {4,3,1}. The complement of this partition is
 * given by ~13, equal to 0010 (integer 2) (after masking leading 1's).
 */

public class Wheatley {
    
    final int n;
    final double[] fQ, fC;
    final double[][] fA, fB;
    
    public Wheatley(int n) {
        this.n = n;
        int nf = 1<<n;  // 2^n
        fQ = new double[nf];
        fC = new double[nf]; // we actually only use half the storage
        for(int i=0; i<n; i++) {
            fQ[1<<i] = 1.0;
        }
        fA = new double[n][nf];
        fB = new double[n][nf];
    }
    
    /*
     * Computation of sum of connected diagrams.
     */
    public double fcCalc(double[][] eArray) {
        int nf = 1<<n;
        
        for (int i=3; i<nf; i++) {
            int j = i & -i;
            if (i==j) continue; // 1-point set
            int k = (i&~j); // k is the points in i other than j
            int jj = Integer.numberOfTrailingZeros(j); // jj = log2(j)
            fQ[i] = fQ[k];
            for (int l=0; l<n; l++) {
                int ll = 1<<l;
                if ((ll&k)==0) continue;
                // l is a point in i, but is not j
                fQ[i] *= eArray[jj][l];
            }
        }

        //Compute the fC's
        for(int i=1; i<nf; i++) {
            fC[i] = fQ[i];
            int iLowBit = i & -i;
            for(int j=1; j<i; j++) {
                if ((j & iLowBit) == 0) continue;
                int jComp = i & ~j;
                if ((jComp | j) != i) continue;
                fC[i] -= fC[j] * fQ[jComp];//for fQ, flip the bits on j; use only those appearing in i
            }
        }

        // find fA1
        for (int i=1; i<nf; i++) {
            fA[0][i] = 0;
            int jBits = (i&-i);
            if (jBits != 1) {
                jBits |= 1;
            }
            else {
                int ii = i & ~jBits;
                if (ii==0) {
                    continue;
                }
                jBits |= (ii & ~ii);
            }
            for (int j=3; j<i; j+=2) {
                if ((j & jBits) != jBits) continue;
                int jComp = i & ~j;
                if (jComp==0 || (jComp | j) != i) continue;
                fA[0][i] += fB[0][j] * fC[jComp|1];
            }
            fB[0][i] = fC[i] - fA[0][i];
        }
        
        for (int v=1; v<n; v++) {
            int vs1 = 1<<v;
            for (int i=1; i<nf; i++) {
//                if ((i & vs1) == 0) continue;
                fA[v][i] = 0;
                int jBits = (i&-i);
                if (jBits != vs1) {
                    jBits |= vs1;
                }
                else {
                    int ii = i & ~jBits;
                    if (ii==0) {
                        fB[v][i] = fB[v-1][i] - fA[v][i];
                        continue;
                    }
                    jBits |= (ii & ~ii);
                }
                for (int j=3; j<i; j++) {
                    if ((j & jBits) != jBits) continue;
                    int jComp = i & ~j;
                    if (jComp==0 || (jComp | j) != i) continue;
                    fA[v][i] += fB[v][j] * (fB[v][jComp|vs1] + fA[v][jComp|vs1]);
                }
                fB[v][i] = fB[v-1][i] - fA[v][i];
            }
        }

        return fB[n-1][nf-1];
    }
    
    public static void main(String[] args) {
        int n = 4;
        Wheatley w = new Wheatley(n);
        double[][] eArray = new double[n][n];
        for(int k1=0; k1<n; k1++) {
            for(int k2=0; k2<k1; k2++) {
                eArray[k1][k2] = (k1+2)+(k2+2);//just fill in some values for development
                eArray[k2][k1] = eArray[k1][k2];
                System.out.println(k1+" "+k2+" "+eArray[k1][k2]);
            }
        }
        double fc = w.fcCalc(eArray);
        System.out.println("fc "+fc);
    }

}
