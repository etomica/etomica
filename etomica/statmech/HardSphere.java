package etomica.statmech;

/**
 * Hard sphere properties based on the Hall EOS.
 */
public class HardSphere {
    
    private double sigma = 1.0;
    
    public HardSphere(){
    }
            
    public double packingFraction(double rho){
        return Math.PI*(sigma*sigma*sigma)*rho/6.0;
    }
    
    public double zSolid(double y){
        double a = y*y;//y^2
        double b = a*y;//y^3
        double c = b*y;//y^4
        double d = c*y;//y^5
        double e = d*y;//y^6
        double numerator1  = (1.0 + y + a - 0.67825*b - c -0.5*d);
        double numerator2  = (6.028*e*Math.exp(((Math.PI*Math.sqrt(2.0)/6.0) - y)*(7.9 - 3.9*((Math.PI*Math.sqrt(2.0)/6.0) - y))));
        double denominator = (1.0 - 3.0*y + 3.0*a - 1.04305*b);
        return (numerator1 - numerator2)/denominator;
    }
    
    
    public double integrand(double rho){
        double y = packingFraction(rho);
        double z = zSolid(y);
        return (z - 1.0)/rho;
    }
        
    
    
    public double doNumericalIntegration(double a,double b,int n){
        double h = (b - a)/(double)n;
        int k = 1;
        double sum = 0.0;
        while(k <= n){
            sum += integrand(a+k*h);
            k++;
        }
        return (h*sum);    
    }
    
    public double idFreeEnergy(double rho){return (Math.log(rho)-1.0);}
    
    public void setSigma(double s){sigma = s;}
    
    public static void main(String[] args) {

        HardSphere hs = new HardSphere();
        double density = 1.3;
        double fid = hs.idFreeEnergy(1.04086);
        double fex = hs.doNumericalIntegration(1.04086,density,100000);
        System.out.println("Ideal Gas Free Energy = "+fid);
        System.out.println("Excess Free Energy = "+fex);
        System.out.println("Absolute Free Energy = "+(fex+fid+5.8644));//Fexcess = 5.8644 for N = 32;5.9117 for N = 108 ;5.9208 for N = 256 JCP81 Pg 3191
    }
}