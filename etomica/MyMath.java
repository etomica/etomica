package simulate;

public class MyMath {
    
    public static int pow(int k, int n) {  //returns k^n (k raised to the nth power), with 0^0 defined as 1
        int result = 1;
        for(int i=n; --i>=0;) {result *= k;}
        return result;
    }
    
    public static class Array {
        
        /**
         * Generate a vector of n integers all having value x.
         */
        public static int[] fillInt(int n, int x) {
            int[] a = new int[n];
            for(int i=0; i<n; i++) {a[i] = x;}
            return a;
        }
        
        /**
         * Generate a n-by-n array with all diagonal elements equal to x (off-diagonal elements zero)
         */
        public static double[][] diagonalMatrix(int n, double x) {
            double[][] a = new double[n][n];
            for(int i=0; i<n; i++) {a[i][i] = x;}
            return a;
        }
    }
}