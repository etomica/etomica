package etomica.virial;

/* 
 * Implementation of Wheatley's algorithm for computing sum of biconnected diagrams. For a graph of order n
 * the partitions are represented by bit strings of length n, formed as integers from 1 to 2^n. Thus, for example,
 * for n = 4, the bits 1101 (integer 13) indicates the partition {4,3,1}. The complement of this partition is
 * given by ~13, equal to 0010 (integer 2) (after masking leading 1's).
 */

public class Wheatley {
    
    final int n, nf, mask;
    final int[] iMask;
    final double[] fQ, fC;
    
    public Wheatley(int n) {
        this.n = n;
        nf = 1<<n;  // 2^n
        mask = 2*nf - 1;//not used
        System.out.println(nf);
        fQ = new double[nf];
        fC = new double[nf];
        iMask = new int[n];
        for(int i=0; i<n; i++) {
            iMask[i] = 1<<i;  //all bits 0 except for bit i
            fQ[iMask[i]] = 1.0;
 //           fC[iMask[i]-1] = 1.0;
        }
    }
    
    /*
     * Computation of sum of connected diagrams.
     */
    public double fcCalc(double[][] eArray) {
        
        //Compute all of the fQ's
        for(int k1=1; k1<n; k1++) {
            for(int k2=0; k2<k1; k2++) {
                int index = (iMask[k1] | iMask[k2]);
                System.out.println(index);
                fQ[index] = eArray[k1][k2];
            }
        }
        
        //Fill fQ with values of prod_j[e_{Mj}], where M is most significant bit in index and
        //product is over all j having nonzero bit in index
        for(int k1=1; k1<n; k1++) {
            for(int k2=0; k2<k1; k2++) {
                double e = eArray[k1][k2];
                for(int i=0; i<iMask[k2]; i++) {//loop over indices below where k2 is most significant bit
                    int index1 = i | iMask[k1];
                    int index2 = index1 | iMask[k2];
                    System.out.println(index2);
                    fQ[index2] = e * fQ[index1];
                }
            }
        }
        
        //Fill fQ with values of prod_{ij}[e_{ij}], where product is over all i and j
        //having nonzero bit in index. This is done by taking prod_[Mj} times prod_{ij}
        //where i < M
        for(int k1=1; k1<n; k1++) {
            for(int i=1; i<iMask[k1]; i++) {
                fQ[i | iMask[k1]] *= fQ[i];
            }
        }
        
        //Compute the fC's
        for(int i=1; i<nf; i++) {
            fC[i] = fQ[i];
            int iLowOneBit = i & -i;
            for(int j=1; j<i; j++) {
                //if((i & j) == 0) continue; //see that i and j have some bits in common
                if ((iLowOneBit & j) == 0) continue;
                int jComp = i & ~j;
                if ((jComp | j) != i) continue;
                fC[i] -= fC[j] * fQ[jComp];//for fQ, flip the bits on j; use only those appearing in i
                System.out.println(i+" "+fC[i]);
            }
        }
        
        return fC[nf-1];
    }
    
    public static void main(String[] args) {
        int n = 3;
        Wheatley w = new Wheatley(n);
        double[][] eArray = new double[n][n];
        for(int k1=0; k1<n; k1++) {
            for(int k2=0; k2<k1; k2++) {
                eArray[k1][k2] = (k1+2)+(k2+2);//just fill in some values for development
            }
        }
        double fc = w.fcCalc(eArray);
        System.out.println("fc "+fc);
        int maxInt = (1<<31) - 1;
        int mask = (1<<4)-1;
        System.out.println(mask & (~13));
        System.out.println(maxInt);
    }

}
