package etomica.virial;

import etomica.potential.P2LennardJones;
import etomica.space.Space;
import etomica.space3d.Space3D;

/**
 * 
 * This is just a main method for PercusYevick.java...
 * 
 * It creates a discretization of the Lennard-Jones Mayer function, which PercusYevick then employs to compute Percus-Yevick
 * approximations of the virial coefficients up to (m+1) order.  The second and third coefficients are fully accurate (the 
 * PY approximation is exact).
 * 
 * @author kate
 *
 */


public class BnPYLJ {
	
public static void main(String[] args) {
        
		System.out.println("Literature values for sigma = 1 and T* = 1");
		System.out.println("B2 = -5.3158 (Sun & Teja 1996)");
		System.out.println("B3 =  1.8849 (Sun & Teja 1996) ");
		System.out.println("B4PY = -2.9394 (Dyer et al 2001)");
		System.out.println("B5PY = -37.402 (Barker et al 1966)\n");
		
		int power = 20; // Defines discretization
		double reducedTemp = 1.0; // kT/epsilon
		double r_max = 50; // Defines range of separation distance, r = [0 rmax]
		
		if (args.length == 0) {
		}
		else if (args.length == 2) {
			power = Integer.parseInt(args[0]);
			r_max = Double.parseDouble(args[1]);
		} else {
        	throw new IllegalArgumentException("Incorrect number of arguments passed.");
        }
		
		// Number of grid points in the discretizations of r- and k-space
		int N = (int) Math.pow(2, power) - 1;
    	
		// m+1 is the highest order of virial coefficient to be calculated
        int m = 5; 
        
        double del_r = r_max/(N-1);

		double r = 0.0;
	
		double[] fr = getfr( N, del_r,reducedTemp);
		
		PercusYevick py = new PercusYevick(); 
		double[] B = py.computeB(fr, m, N, del_r);
		
		System.out.println("Values computed here at T* = " + reducedTemp + ":");
		
		for (int i=0;i<m;i++) {
			System.out.println("B" + (i+2) + " = " + B[i]);
		}
		

	}

	public static double[] getfr(int N, double del_r, double reducedTemp) {
		
		// Lennard-Jones Potential
	    double sigma = 1.0;
		double epsilon = 1.0;
		Space space = Space3D.getInstance();
		P2LennardJones p2 = new P2LennardJones(space, sigma, epsilon);
		double u;
		
		double r = 0.0;
		
		double[] fr = new double[N];  // Holds discretization of Mayer function in r-space
		
		for (int n = 0; n<N; n++) {
			
			u = p2.u(r*r);
			
			double x = -u/reducedTemp;
			
			if ( Math.abs(x) < 0.01) {
				
				fr[n] = x + x*x/2.0 + x*x*x/6.0 + x*x*x*x/24.0 + x*x*x*x*x/120.0;
				
				
			} else {
				
				fr[n] = Math.exp(x)-1.0;
				
			}
			
			// fr[n] = Math.exp(x)-1.0;
			
			
			r += del_r; 
	
		}
		
		return fr;
		
	}


}
