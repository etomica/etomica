package etomica.virial;

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
            fQ[iMask[i]-1] = 1.0;
 //           fC[iMask[i]-1] = 1.0;
        }
    }
    
    public double fcCalc(double[][] eArray) {
        
        //Compute all of the fQ's
        for(int k1=1; k1<n; k1++) {
            for(int k2=0; k2<k1; k2++) {
                int index = (iMask[k1] | iMask[k2]);
                System.out.println(index-1);
                fQ[index-1] = eArray[k1][k2];
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
                    System.out.println(index2-1);
                    fQ[index2-1] = e * fQ[index1-1];
                }
            }
        }
        
        //Fill fQ with values of prod_{ij}[e_{ij}], where product is over all i and j
        //having nonzero bit in index. This is done by taking prod_[Mj} times prod_{ij}
        //where i < M
        for(int k1=1; k1<n; k1++) {
            for(int i=1; i<iMask[k1]; i++) {
                fQ[(i | iMask[k1]) - 1] *= fQ[i-1];
            }
        }
        
        //Compute the fC's
        for(int i=0; i<nf; i++) {
            fC[i] = fQ[i];
            for(int j=0; j<i; j++) {
                if((i & j) == 0) continue; //see that i and j have some bits in common
                fC[i-1] -= fC[j-1] * fQ[(i & ~j)-1];//flip the bits on j; use only those appearing in i
            }
        }
        
        return fC[nf-1];
    }
    
    public static void main(String[] args) {
        Wheatley w = new Wheatley(5);
        double[][] eArray = new double[5][5];
        for(int k1=0; k1<5; k1++) {
            for(int k2=0; k2<k1; k2++) {
                eArray[k1][k2] = k1+k2;//just fill in some values for development
            }
        }
        double fc = w.fcCalc(eArray);
    }

}
