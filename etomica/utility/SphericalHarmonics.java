package etomica.utility;


public class SphericalHarmonics{
 
    private AssociatedPolynomial ap ;
    
    public SphericalHarmonics(){
        ap = new AssociatedPolynomial();
    }
        
    private double factor(int l,int m){
        
        return Math.sqrt(((2.0*l + 1.0)*SpecialFunctions.factorial(l-m))/(4.0*Math.PI*SpecialFunctions.factorial(l+m)));
        
    }
    
    private double y(int l,int m,double theta){
        return factor(l,m)*ap.plgndr(l, m, theta);
    }
    
    public double realYm(int l,int m,double theta,double phi){
        if(m < 0.0){
            m = Math.abs(m);
            return Math.pow(-1,m)*y(l,m,theta)*Math.cos(m*phi);
        }
        else{
            return y(l,m,theta)*Math.cos(m*phi);
        }
    }
    
    public double imaginaryYm(int l,int m,double theta,double phi){
       if(m < 0.0){
            m = Math.abs(m);
            return Math.pow(-1,m)*y(l,m,theta)*Math.sin(m*phi)*(-1.0);
       }
       else{ 
         return y(l,m,theta)*Math.sin(m*phi);
       }
    }
    
    
}