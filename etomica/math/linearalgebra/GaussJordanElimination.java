package etomica.math.linearalgebra;

//maybe wrong
//who wrote this?
public class GaussJordanElimination implements java.io.Serializable {
    
    private GaussJordanElimination(){}
    
    /**
     * Performs Gauss-Jordan elimination in place to solve Ax=b.  A is
     * an nxn matrix, b is an nxm vector.
     */
    public static void gaussj(double [][]a,int n,double [][]b,int m){
        
        int[] indxc = new int[n];
        int[] indxr = new int[n];
        boolean[] ipiv  = new boolean[n];
        
        for(int j = 0; j < n; j++)
            ipiv[j] = false;
        
        for(int i = 0 ; i < n; i++){
            double max = -1;
            int icol = -1;
            int irow = -1;
            
            // loop over rows
            for(int j = 0; j < n; j++) {
                // if we haven't already pivotted this row
                if(!ipiv[j]) {
                    // loop over cols
                    for(int k = 0; k < n; k++){
                        // if we haven't already pivotted the row we'd swap with
                        // if this was the max
                        if(!ipiv[k]){
                            if(Math.abs(a[j][k]) >= max){
                                max = Math.abs(a[j][k]);
                                irow = j;
                                icol = k;
                            }
                        }
                    }
                }
            }

            if (max == 0) {
                // singular matrix.
                for (int j=0; j<n; j++) {
                    for (int k=0; k<m; k++) {
                        b[j][k] = Double.NaN;
                    }
                }
                return;
            }
                    
            // mark the column as having been pivotted
            ipiv[icol] = true;

            if(irow != icol){
                //swap rows irow and icol
                for(int l = 0;l < n; l++) {
                    double temp = a[irow][l];
                    a[irow][l] = a[icol][l];
                    a[icol][l] = temp;
                }
                
                for(int l = 0;l < m; l++) {
                    double temp = b[irow][l];
                    b[irow][l] = b[icol][l];
                    b[icol][l] = temp;
                }
            }

            // remember pivotting scheme... for the ith step, we swapped
            // irow and icol!
            indxr[i] = irow;
            indxc[i] = icol;

            double pivinv = 1.0/a[icol][icol];
            // divide the icol row by the pivot value
            for(int l = 0;l < n;l++) {
                a[icol][l]*= pivinv;
            }
            for(int l = 0;l < m;l++) {
                b[icol][l]*= pivinv;
            }

            // subtract the icol row from the other rows (times a constant so 
            // that the icol element of a becomes 0).
            for(int ll = 0;ll < n;ll++) {
                if(ll != icol){
                    double dum = a[ll][icol];
                    for(int l = 0;l < n; l++){
                        a[ll][l] -= a[icol][l]*dum;
                    }
                    for(int l = 0;l < m; l++){
                        b[ll][l] -= b[icol][l]*dum;
                    }  
                }
            }
        }
 
        for(int l = n-1;l >= 0;l--){
            // unswap the rows back to where they came from
            if(indxr[l] != indxc[l]) {
                int irl = indxr[l];
                int icl = indxc[l];
                for(int k = 0;k < n;k++) {
                    double temp = a[k][irl];
                    a[k][irl] = a[k][icl];
                    a[k][icl] = temp;
                }
            }
        }
    }
    
    public static void main(String[] args) {
        double[][] aa = new double[3][3];
        double[][] bb = new double[3][1];
        aa[0][0] = 2.0;
        aa[0][1] = 1.0;
        aa[0][2] = 1.0;
        aa[1][0] = 4.0;
        aa[1][1] = 1.0;
        aa[1][2] = 0.0;
        aa[2][0] = -2.0;
        aa[2][1] = 2.0;
        aa[2][2] = 1.0;
        bb[0][0] = 1.0;
        bb[1][0] = -2.0;
        bb[2][0] = 7.0;
        
        GaussJordanElimination.gaussj(aa, 3, bb, 1);

        for(int i = 0;i < 3;i++){
            System.out.println(bb[i][0]+"\t"+i);
        }
    }
    
}
