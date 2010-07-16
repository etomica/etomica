package etomica.virial;

import etomica.potential.P2LennardJones;
import etomica.space.Space;
import etomica.space3d.Space3D;

/**
 * 
 * This is just a main method for PercusYevick.java...
 * 
 * It creates a discretization of the Lennard-Jones Mayer function, which PercusYevick then employs to compute Percus-Yevick
 * approximations of the virial coefficients up to mth order.  The second and third coefficients are fully accurate (the 
 * PY approximation is exact).
 * 
 * @author kate
 *
 */


public class BnPYLJ {
	
public static void main(String[] args) {
        
		// To make sure that everything is working fine:
		System.out.println("Literature values for LJ with sigma = 1 and T* = 1");
		System.out.println("B2 = -5.3158 (Sun & Teja 1996)");
		System.out.println("B3 =  1.8849 (Sun & Teja 1996) ");
		System.out.println("B4PY = -2.9394 (Dyer et al 2001)");
		System.out.println("B5PY = -37.402 (Barker et al 1966)\n");
		
		
		int power = 20; // Defines discretization
		double reducedTemp = 1.0; // kT/epsilon
		double r_max = 150; // Defines range of separation distance, r = [0 rmax]
		
		if (args.length == 0) {
		}
		else if (args.length == 2) {
			power = Integer.parseInt(args[0]);
			r_max = Double.parseDouble(args[1]);
		} else {
        	throw new IllegalArgumentException("Incorrect number of arguments passed.");
        }
		
		// Number of grid points in the discretizations of r- and k-space
		int N = (int) Math.pow(2, power);
		double del_r = r_max/(N-1);
		
		// Get Mayer function for this discretization
		double[] fr = getfr( N, del_r,reducedTemp);
	
		
        int m = 6; // highest order of virial coefficient to be calculated
        
		PercusYevick py = new PercusYevick(); 
		double[] B = py.computeB(fr, m, N, del_r, false);
		
		System.out.println("Values computed here at T* = " + reducedTemp + ":");
		
		for (int i=2;i<=m;i++) {
			System.out.println("B" + i + " = " + B[i-2]);
		}
		

	}

	public static double[] getfr(int N, double del_r, double reducedTemp) {
	    
		double sigma = 1.0;
		
	    double epsilon = 1.0;
		
		Space space = Space3D.getInstance();
		
		P2LennardJones p2 = new P2LennardJones(space, sigma, epsilon);
		
		double r = 0.0;
		
		double[] fr = new double[N];  // Holds discretization of Mayer function in r-space
		
		for (int n = 0; n<N; n++) {
			
			double u = p2.u(r*r);
			
			double x = -u/reducedTemp;
			
			if ( Math.abs(x) < 0.01) {
				
				fr[n] = x + x*x/2.0 + x*x*x/6.0 + x*x*x*x/24.0 + x*x*x*x*x/120.0;	
				
			} else {
				
				fr[n] = Math.exp(x)-1.0;
				
			}
			
			r += del_r; 
	
		}
		
		return fr;
		
	}


}
