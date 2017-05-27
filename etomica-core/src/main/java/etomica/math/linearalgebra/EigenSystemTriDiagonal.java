/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.math.linearalgebra;


//Eigenvalues and Eigenvectors of a tridiagonal matrix Refer Numerical Recipes page 480-481 
public class EigenSystemTriDiagonal implements java.io.Serializable {
    
    private int n = 20;
    private double []d ;//contains diagonal elements of the tridiagonal matrix
    private double []e ;//subdiagonal elements of the tridiagonal matrix with e[1] arbitrary
                                        //on output e is destroyed
    private double [][]z ;//identity matrix                                    
    
        
    public EigenSystemTriDiagonal(){
    }
    
    
    public void tqli(double []diagonal,double []offdiagonal,int N){
        n = N;
        d = new double [n];
        e = new double [n];
        z = new double [n][n];
       for(int i = 0;i < n;i++){ 
            d[i] = diagonal[i];
            e[i] = offdiagonal[i];
       }
       for(int i = 0;i < n;i++){
        for(int j = 0;j < n;j++){
            if(i != j) {z[i][j] = 0.0;}
            else{
                //i == j
                z[i][j] = 1.0;
            }
        }
       }
        int     m,l,iter,i,k;
        double  s,r,p,g,f,dd,c,b;
        
        for(i = 1;i < n;i++){
            e[i-1] = e[i];
        }
            e[n-1] = 0.0;
        for(l = 0;l  < n;l++){
            iter = 0;
            do {
                for(m = l;m <(n-1);m++){
                    dd = Math.abs(d[m]) + Math.abs(d[m+1]);
                    if(Math.abs(e[m] + dd) == dd)break;
                }
                if(m != l){
                    if(iter++ == 30){
                        System.err.println("Too many iterations in tqli");
                    }
                    g = (d[l+1] - d[l])/(2.0*e[l]);
                    r = pythag(g,1.0);
                    g = d[m] - d[l] + e[l]/(g + sign(r,g));
                    s = c = 1.0;
                    p = 0.0;
            for(i = m-1;i >= l;i--){
                f = s*e[i];
                b = c*e[i];
                e[i+1] = ( r = pythag(f,g));
                if ( r == 0.0){
                    d[i+1] -= p;
                    e[m]    = 0.0;
                    break;
                }
                s = f/r;
                c = g/r;
                g = d[i+1] - p;
                r = (d[i] - g)*s + 2.0*c*b;
                d[i+1] = g + (p = s*r);
                g = c*r - b;
                //Eigen Vector Calculation
                for(k = 0;k < n;k++){
                    f = z[k][i+1];
                    z[k][i+1] = s*z[k][i] + c*f;
                    z[k][i]   = c*z[k][i] - s*f;//System.out.println(z[k][i]);//eigen vectors
                }
            }
            if(r == 0.0 && i >= 1) continue;
            d[l] -= p;
            e[l]  = g;
            e[m]  = 0.0;
                }
            }while( m != l);System.out.println("Eigen value:"+l+"\t"+d[l]);//eigen values
        }
           output(); 
    }
    
    public double getEigenValue(int i){
        return d[i];
    }
    
    public double getEigenVectorElement(int i, int j){
        return z[i][j];
    }
    
    public double [][]getEigenVectors(){
        double [][]r = new double [n][n];
        for(int i = 0; i < n;i++){
            for(int j = 0;j < n;j++){
                r[i][j] = z[i][j];
            }
        }
        return r;
    }
        
    
    public void output(){
        for(int i = 0; i < n;i++){
            for(int j = 0;j < n;j++){
                System.out.println("Eigen vector:"+"\t"+i+"\t"+j+"\t"+z[i][j]);
            }
        }
    }
    
    public double sign(double a , double b){
        if(b > 0.0){
            return Math.abs(a);
        }
        return -Math.abs(a);
    }
    
    public double pythag(double a,double b){
        double absa , absb;
        absa = Math.abs(a);
        absb = Math.abs(b);
        if(absa > absb){
            return absa*Math.sqrt(1.0 + sqr(absb/absa));
        }
        return (absb == 0.0 ? 0.0:absb*Math.sqrt(1.0 + sqr(absa/absb)));
    }
    
    public double sqr(double a){
        return a*a;
    }
}
