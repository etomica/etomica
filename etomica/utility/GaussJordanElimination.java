package etomica.utility;


public class GaussJordanElimination{
    
    private int []indxc;
    private int []indxr;
    private int []ipiv;
    
    private int    i,icol,irow,j,k,l,ll;
    private double big,dum,pivinv,temp;
    
    public double [][]aa = new double [3][3];
    public double [][]bb = new double[3][1];
    
    public GaussJordanElimination(){}
    
    
    public void gaussj(double [][]a,int n,double [][]b,int m){
        
        indxc = new int [n];
        indxr = new int [n];
        ipiv  = new int [n];
        
        for(j = 0 ;j < n;j++)
            ipiv[j] = 0;
        
        
        for(i = 0 ; i < n; i++){
            big = 0.0;
            for(j = 0;j < n;j++)
                if(ipiv[j] != 1.0)
                for(k = 0;k < n;k++){
                    if(ipiv[k] == 0){
                        if(Math.abs(a[j][k]) >= big){
                            big = Math.abs(a[j][k]);
                            irow = j;
                            icol = k;
                        }
                    }
                    else{
                        if(ipiv[k] > 1){System.err.println("gaussj:Singular Matrix-1");}
                    }
                
            }
        
                    ++(ipiv[icol]);
                    
                    if(irow != icol){
                        for(l = 0;l < n; l++)
                            sWAP(a[irow][l],a[icol][l]);
                        
                        for(l = 0;l < m; l++)
                            sWAP(b[irow][l],b[icol][l]);
                        
                    }
                    
                    indxr[i] = irow;
                    indxc[i] = icol;
                    if(a[icol][icol] == 0.0){System.err.println("gaussj:Singular Matrix-2");}
                    pivinv = 1.0/a[icol][icol];
                    a[icol][icol] = 1.0;
                    for(l = 0;l < n;l++)a[icol][l]*= pivinv;
                    for(l = 0;l < m;l++)b[icol][l]*= pivinv;
                    
                    for(ll = 0;ll < n;ll++)
                        if(ll != icol){
                            dum = a[ll][icol];
                            a[ll][icol] = 0.0;
                            for(l = 0;l < n;l++){a[ll][l] -= a[icol][l]*dum;}
                            for(l = 0;l < m;l++){b[ll][l] -= b[icol][l]*dum;}  
                        }
                    }
         
                    for(l = n-1;l >= 0;l--){
                        if(indxr[l] != indxc[l])
                            for(k = 0;k < n;k++)
                                sWAP(a[k][indxr[l]],a[k][indxc[l]]);
                            
                        
                    }
       }
    
    public void sWAP(double a, double b){
        double temp = a;
        a = b;
        b = temp;
    }
    
    public void definition(){
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
    }
    
        public void output(){
            for(int i = 0;i < 3;i++){
             System.out.println(bb[i][0]+"\t"+i);
            }
        }
    
    public static void main(String[] args) {
        GaussJordanElimination gje = new GaussJordanElimination();
        gje.definition();
        gje.gaussj(gje.aa,3,gje.bb,1);
        gje.output();
    }
    
}
