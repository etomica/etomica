package etomica.utility;


public class SphericalHarmonics{
             
    private static double factor(int l,int m){
        return Math.sqrt(((2.0*l + 1.0)*SpecialFunctions.factorial(l-m))/(4.0*Math.PI*SpecialFunctions.factorial(l+m)));
    }
    
    private static double y(int l,int m,double theta){
        return factor(l,m)*AssociatedPolynomial.plgndr(l, m, theta);
    }
    
    public static double realYm(int l, int m, double theta, double phi){
        if(m < 0.0){
            m = -m;
            return Math.pow(-1,m)*y(l,m,theta)*Math.cos(m*phi);
        }
        else {
            return y(l,m,theta)*Math.cos(m*phi);
        }
    }
    
    public static double imaginaryYm(int l,int m,double theta,double phi){
       if(m < 0.0){
            m = -m;
            return Math.pow(-1,m)*y(l,m,theta)*Math.sin(m*phi)*(-1.0);
       }
       else { 
         return y(l,m,theta)*Math.sin(m*phi);
       }
    }
}