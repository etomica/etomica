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
            int j = i & -i;//lowest bit in i
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
            fB[0][i] = fC[i];
            if((i & 1) == 0 || i == 1) continue;//if i doesn't contain 1, or is 1, fA and fB are done
            //at this point we know lowest bit in i is 1; add next lowest bit to it
            int ii = i - 1;//all bits in i but lowest
            int jBits = 1 | (ii & -ii);
            //at this point jBits has 1 and next lowest bit in i
            for (int j=3; j<i; j+=2) {//sum over partitions of i containing 1
                if ((j & jBits) != jBits) continue;//ensure jBits are in j
                int jComp = i & ~j; //subset of i complementing j
                if (jComp==0 || (jComp | j) != i) continue;
                fA[0][i] += fB[0][j] * fC[jComp|1];
            }
            fB[0][i] -= fA[0][i];//remove from B graphs that contain articulation point at 0
        }
        
        for (int v=1; v<n; v++) {
            int vs1 = 1<<v;
            for (int i=1; i<nf; i++) {
                fA[v][i] = 0;
                fB[v][i] = fB[v-1][i];
                if ((i & vs1) == 0) continue;//if i doesn't contain v, fA and fB are done
                int jBits = (i&-i); //lowest bit in i
                if (jBits != vs1) { //lowest bit is not v; add v to it
                    jBits |= vs1;
                }
                else { //lowest bit is v; add next lowest bit to it
                    int ii = i & ~jBits;
                    if (ii==0) { //lowest bit is only bit
                        continue;
                    }
                    jBits |= (ii & -ii);
                }
                //at this point jBits has (lowest bit + v) or (v + next lowest bit)
                for (int j=3; j<i; j++) {//sum over partitions of i
                    if ((j & jBits) != jBits) continue;//ensure jBits are in j
                    int jComp = i & ~j;//subset of i complementing j
                    if (jComp==0 || (jComp | j) != i) continue;
                    fA[v][i] += fB[v][j] * (fB[v][jComp|vs1] + fA[v][jComp|vs1]);
                }
                fB[v][i] -= fA[v][i];//remove from B graphs that contain articulation point at v
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
