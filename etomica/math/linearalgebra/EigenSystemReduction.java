package etomica.math.linearalgebra;

//Reduction of a symmetric matrix to tridiagonal form Refer Numerical Recipes page474-475 
public class EigenSystemReduction {
    
    private int n = 20;
    private double []d ;//contains diagonal elements of the tridiagonal matrix
    private double []e ;//subdiagonal elements of the tridiagonal matrix with e[1] arbitrary
                                        //on output e is destroyed
    private double [][]a ;//input symmetric matrix to be reduced to tridiagonal form                                    
    
        
    public EigenSystemReduction(){
    }
    
    
    public void tred2(double [][]A,int N){
        n = N;
        a = new double [n][n];
        d = new double [n];
        e = new double [n];
        for(int m = 0; m < N;m++){
            for(int l = 0;l < N;l++){
                a[m][l] = A[m][l];
            }
        }
        
        int l,k,j,i;
        double scale,hh,h,g,f;
        for(i = n-1;i >= 1;i--){
            l = i-1;
            h = scale = 0.0;
            if(l > 1){
                for(k = 0; k <= l;k++)
                    scale += Math.abs(a[i][k]);
                if(scale == 0.0)
                    e[i] = a[i][l];
                else{
                    for(k = 0;k <= l;k++){
                        a[i][k] /= scale;
                        h += a[i][k]*a[i][k];
                    }
                    f = a[i][l];
                    g = (f >= 0.0 ? -1*Math.sqrt(h):Math.sqrt(h));
                    e[i] = scale*g;
                    h -= f*g;
                    a[i][l] = f - g;
                    f = 0.0;
                    for(j = 0;j <= l;j++){
                        a[j][i] = a[i][j]/h;
                        g = 0.0;
                        for(k = 0;k <= l;k++)
                            g += a[j][k]*a[i][k];
                        for(k = j+1;k <= l;k++)
                            g += a[k][j]*a[i][k];
                        e[j] = g/h;
                        f += e[j]*a[i][j];
                    }
                    hh = f/(h+h);
                    for(j = 0;j <= l;j++){
                        f = a[i][j];
                        e[j] = g = e[j] - hh*f;
                        for(k = 1;k <= j;k++)
                            a[j][k] -= (f*e[k] + g*a[i][k]);
                    }
                  }
                }else
                    e[i] = a[i][l];
                d[i] = h;
        }
        d[0] = 0.0;
        e[0] = 0.0;
        
        for(i = 0;i < n;i++){
            l = i;
            if(d[i] != 0.0){
                for(j = 0;j <= l;j++){
                    g = 0.0;
                    for(k = 0;k <= l;k++)
                        g += a[i][k]*a[k][j];
                    for(k = 0;k <= l;k++)    
                        a[k][j] -= g*a[k][i];
                }
            }
            d[i] = a[i][i];
            a[i][i] = 1.0;
            for(j = 0;j <= l;j++)
                a[j][i] = a[i][j] = 0.0;
        }output();
    }
                    
                        
    public double []getDiagonalElements(){
                  return d;
             }
        
    public double []getOffDiagonalElements(){
                  return e;
             }
    
    
    public void output(){
        for(int i = 0; i < n;i++){
           System.out.println(d[i]+"\t"+e[i]);
        }
    }
    
    
    
}