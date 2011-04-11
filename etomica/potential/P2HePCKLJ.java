
package etomica.potential;
import etomica.space.ISpace;
import etomica.space.Space;
import etomica.space3d.Space3D;
import etomica.units.Dimension;
import etomica.units.Energy;
import etomica.units.Length;

/**
 * Ab initio pair potential for helium from Przybytek et al. (2010) Phys. Rev. Lett. 104, 183003.  This is a true pair potential, rather than a pairwise-additive potential.
 * 
 * In this class, only the pair potential is valid, not the gradients, etc.  I am unlikely to ever include those...
 *
 * @author Kate Shaul
 */
public class P2HePCKLJ extends Potential2SoftSpherical {
    
    public P2HePCKLJ(ISpace space) {
        super(space);
   
    }

    /**
     * The energy u.
     */
    public double u(double r2) {
    	
    	double r = Math.sqrt(r2);
        
    	r = r/AngstromPerBohrRadius; // Bohr radius
    	
    	
    	//Parameters
    	
    	double[] C = new double[17];
    	C[3] = 0.000000577235; 
    	C[4] = -0.000035322; 
    	C[5] = 0.000001377841; 
    	C[6] = 1.461830; 
    	C[8] = 14.12350; 
    	C[10] = 183.7497; 
    	C[11] = -0.7674*100; 
    	C[12] = 0.3372*1e4; 
    	C[13] = -0.3806*1e4; 
    	C[14] = 0.8534*1e5; 
    	C[15] = -0.1707*1e6;  
    	C[16] = 0.286*1e7;  

    	double a = 3.64890303652830;
    	double b = 2.36824871743591;
    	double eta = 4.09423805117871;
    	
    	double[] P = new double[] { -25.4701669416621, 269.244425630616, -56.3879970402079};
    	double[] Q = new double[] {38.7957487310071, -2.76577136772754};

    	

    	double u1 = 0;
    	for (int i=0; i<=2; i++) {
    		
    		u1 = u1 + P[i]*Math.pow(r,i);
    		
    	}
    	u1 = u1*Math.exp(-a*r);
    	
    	double u2 = 0;
    	for (int i=0; i<=1; i++) {
    		
    		u2 = u2 + Q[i]*Math.pow(r,i);
    		
    	}
    	u2 = u2*Math.exp(-b*r);
    	
    	double u3 = 0;
    	double x = eta*r;
    	double d = 1.0 + x + x*x/2.0;
    	double factorialn = 2.0;
    	for (int n=3; n<=16; n++) {
    		
    		factorialn = factorialn*n;
    		d = d + Math.pow(x, n)/factorialn;
    		
    		double f = 1.0-(Math.exp(-x)*d);
    		
    		u3 = u3 + f*C[n]/Math.pow(r,n);
    		
    	}

    	
 
        double u = u1 + u2 - u3; // Hartree
        
        u = u*KPerHartree; // Kelvin
       
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
    	P2HePCKLJ p2 = new P2HePCKLJ(space);

   
    	//Separations printed in paper:
    	double[] rA = new double[] {3.0,4.0,5.0,5.6,6.0,7.0,12.0}; //Bohr radii
    	double u;
    	
    	
     	for (int i=0;i<7;i++) {
    		double r = rA[i]*AngstromPerBohrRadius; //Angstrom
    		u = p2.u(r*r); // Kelvin
    		
    	
    		System.out.println(r+"  "+u);
    	}
    	
    	/*double r = 1;
    	while (r < 20) {
    		r = r + 0.2; //Angstrom
    		u = p2.u(r*r); // Kelvin
    		double e = Math.exp(-u/1000);
    	
    		System.out.println(r+"  "+u + "  " + e);
    	}*/
    	
    	double r = 7.0*AngstromPerBohrRadius;
    	u = p2.u(r*r)*1000; // milliKelvin
    	System.out.println(r+"  "+3*u );

		
    	
    }
   

    private static final double AngstromPerBohrRadius = 0.529177; // Rounding provided by Pryzbytek et al. 2010
    private static final double KPerHartree = 315774.65; // Rounding provided by Pryzbytek et al. 2010
    private static final long serialVersionUID = 1L;
    
}
