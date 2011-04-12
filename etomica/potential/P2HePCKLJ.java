
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
 * So far, does not include Vret.
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
    	
    	double C6BO = 1.460977837725; //pulled from potentials.f90...
    	
    	double[] A = new double[6];
    	A[1] = 8.454943177941253;
    	A[2] = 15.552891567112597; 
    	A[3] = 7.559160092169168; 
    	A[4] = 1.417737689876350; 
    	A[5] = 0.142506077478301; 
    	
    	double[] B = new double[7];
    	B[1] = 8.454943177941253;
    	B[2] = 16.006586066260556;
    	B[3] = 10.378373954734820;
    	B[4] = 3.515803817223855;
    	B[5] = 0.591502377533792;
    	B[6] = 0.059455768329599;


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
    	double f3 = 1.0;
    	double f4 = 1.0;
    	for (int n=3; n<=16; n++) {
    		
    		factorialn = factorialn*n;
    		d = d + Math.pow(x, n)/factorialn;
    		
    		double f = 1.0-(Math.exp(-x)*d);
    		
    		if (n==3) {
    			f3 = f;
    		} else if (n==4) {
    			f4 = f;
    		}
    		
    		u3 = u3 + f*C[n]/Math.pow(r,n);
    		
    	}

    	
 
        double u = u1 + u2 - u3; // Hartree
        
        double sumA = 0;
        double sumB = 0;
        double alpha = 1.0/137.036; //fsalpha from potentials.f90...
        x = alpha*r;
        double xn=1.0;
        for (int n=1;n<=5;n++) {
        	xn=xn*x;
        	sumA = sumA + A[n]*xn;
        	sumB = sumB + B[n]*xn;
        }
        
        sumB = sumB + B[6]*xn*x;
        
        double g = (1.0+sumA)/(1.0+sumB);
        
        double VCP = -C6BO*Math.pow(r, -6)*g;
        double Vret = VCP + C6BO*Math.pow(r, -6) + 
        C[4]*Math.pow(r, -4) + C[3]*Math.pow(r,-3);
        System.out.println(r+"  "+Vret*KPerHartree);
        
        u = (u+Vret)*KPerHartree; // Kelvin
       
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
    	
    	//Test separations in fortran code:
    	double[] VTotRet = new double [] {286567.288000 ,  104344.430500 ,   36151.809900 ,   11964.292600 ,    3768.978000 ,    1111.080300 ,     292.720500 ,      58.448600 ,      -0.456800 ,      -9.156800 ,     -10.995600 ,      -9.677400 ,      -6.890700 ,      -4.620800 ,      -2.065600 ,      -0.988900 ,      -0.165500 };
    	double[] rTest = new double [] {1.00 , 1.50 , 2.00 , 2.50 , 3.00 , 3.50 , 4.00 , 4.50 , 5.00 , 5.30 , 5.60 , 6.00 , 6.50 , 7.00 , 8.00 , 9.00 , 12.00};
    	double[] VTot = new double [] {286567.283700 ,  104344.428600 ,   36151.808900 ,   11964.292000 ,    3768.977500 ,    1111.080000 ,     292.720300 ,      58.448400 ,      -0.457000 ,      -9.156900 ,     -10.995700 ,      -9.677500 ,      -6.890800 ,      -4.620900 ,      -2.065700 ,      -0.988900 ,      -0.165600};
    	double[] VRet = new double [] {0.004330 ,       0.001890 ,       0.001050 ,       0.000660 ,       0.000450 ,       0.000330 ,       0.000250 ,       0.000190 ,       0.000150 ,       0.000130 ,       0.000120 ,       0.000100 ,       0.000090 ,       0.000070 ,       0.000050 ,       0.000040 ,       0.000020 };
    	
    	double u;
    	
    	System.out.println("VTot retarded (K) ");
    	System.out.println();
    	System.out.println("r(a0) \t test.f90 \t here    ");
     	for (int i=0;i<rTest.length;i++) {
    		double r = rTest[i]*AngstromPerBohrRadius; //Angstrom
    		u = p2.u(r*r); // Kelvin
    		
    	 
    		System.out.println(rTest[i]+"   \t"+VTotRet[i]+"  \t" +u);
    	}
     	
     	System.out.println("VTot (K) ");
    	System.out.println();
    	System.out.println("r(a0) \t test.f90 \t here    ");
     	for (int i=0;i<rTest.length;i++) {
    		double r = rTest[i]*AngstromPerBohrRadius; //Angstrom
    		//p2.setU(1);
    		u = p2.u(r*r); // Kelvin
    		
    	 
    		System.out.println(rTest[i]+"   \t"+VTot[i]+"  \t" +u);
    	}
     	
     	System.out.println("VRet (K) ");
    	System.out.println();
    	System.out.println("r(a0) \t test.f90 \t here    ");
     	for (int i=0;i<rTest.length;i++) {
    		double r = rTest[i]*AngstromPerBohrRadius; //Angstrom
    		u = p2.u(r*r); // Kelvin
    		
    	 
    		System.out.println(rTest[i]+"   \t"+VRet[i]+"  \t" +u);
    	}

    	
    }
   

    private static final double AngstromPerBohrRadius = 0.529177; // Rounding provided by Pryzbytek et al. 2010
    private static final double KPerHartree = 315774.65; // Rounding provided by Pryzbytek et al. 2010
    private static final long serialVersionUID = 1L;
    
}
