package etomica.potential;
import etomica.space.ISpace;
import etomica.space.Space;
import etomica.space3d.Space3D;
import etomica.units.Dimension;
import etomica.units.Energy;
import etomica.units.Length;

/**
 * Pair potential for argon from Aziz (1993) JCP 99(6): 4518.  This is a true pair potential, rather than a pairwise-additive potential.
 * 
 * In this class, only the pair potential is valid, not the gradients, etc.  I am unlikely to ever include those...
 *
 * @author Kate Shaul
 */
public class P2ArgonTangAndToennies2003 extends Potential2SoftSpherical {
    
    public P2ArgonTangAndToennies2003(ISpace space) {
        super(space);
   
    }

    /**
     * The energy u.
     */
    public double u(double r2) {
    	
    	double r = Math.sqrt(r2);
    	r = r/(0.52917720859); // Bohr radius
    	
    	
    	
    	//Parameters
    	double C6 = 64.30;
    	double C8 = 1623;
    	double C10 = 49060;
    	double A =748.3; // Kelvin
    	double b = 2.031; // inverse Angstroms squared
    	
    	double sum;
    	double C2n;
    	
    	double x = b*r;
    	
    	double kfac  = 1;
    	double sum6 = 0;
    	for (int k=0;k<=6;k++) {		
			if (k>0) {kfac=kfac*k;}
			sum6=sum6+Math.pow(x,k)/kfac;
		}
    	
    	double sum8 = sum6;
    	for (int k=7;k<=8;k++) {
			kfac=kfac*k;
			sum8=sum8+Math.pow(x,k)/kfac;
		}
    	

    	double sum10 = sum8;
    	for (int k=9;k<=10;k++) {
			kfac=kfac*k;
			sum10=sum10+Math.pow(x,k)/kfac;
		}
 
    	double f6 = 1-Math.exp(-x)*sum6;
    	double f8 = 1-Math.exp(-x)*sum8;
    	double f10 = 1-Math.exp(-x)*sum10;
    		
 
        double u = A*Math.exp(-b*r) - (f6*C6/Math.pow(r,6)) - (f8*C8/Math.pow(r,8)) - (f10*C10/Math.pow(r,10)); // Kelvin 
        
        u = u*(4.35974417E5/1.3806503); // Kelvin
       
        return u; // Kelvin
    }

    /**
     * The derivative r*du/dr.
     */
    public double du(double r2) {
        
        return 0;
    }

   /**
    * The second derivative of the pair energy, times the square of the
    * separation:  r^2 d^2u/dr^2.
    */
    public double d2u(double r2) {
     
        return 0;
    }
            
    /**
     *  Integral used for corrections to potential truncation.
     */
    public double uInt(double rC) {
        
        return 0;  //complete LRC is obtained by multiplying by N1*N2/V
    }
    
    public static void main(String[] args) {
    	
    	Space space = Space3D.getInstance();
    	P2ArgonTangAndToennies2003 p2 = new P2ArgonTangAndToennies2003(space);
    	
    	//Minimum of the potential:
    	
    	double Rm=7.10*(0.52917720859); //Angstroms
    	double epsilon = 4.54e-4*1e6;  // microHartree (convert to energy through Boltzmann's constant...)
    	
    	System.out.println(Rm+"   "+epsilon);
    	
    	double r = 2.0;
    	double u;
    	double k =  1.3806503e-23; // J/Kelvin
    	double JPerHartree = 4.359743e-18; // J;
    	while (r<10.0) {
    		r = r + 0.05;
    		u = p2.u(r*r); // Kelvin
    		u = u*k/JPerHartree*1e6;// microHartree
    		System.out.println(r+"  "+u);
    	}
    }
   

   
   
    private static final long serialVersionUID = 1L;
    
}
