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
		
		
		int power = 15; // Defines discretization
		
		double r_max = 100; // Defines range of separation distance, r = [0 rmax]
		
		int N = 1<<power;
		double del_r = r_max/(N-1);
		
		int m = 4; // highest order of virial coefficient to be calculated
		double reducedTemp;
		if (args.length == 0) {
		}
		else if (args.length == 4) {
			power = Integer.parseInt(args[0]);
			r_max = Double.parseDouble(args[1]);
			reducedTemp = Double.parseDouble(args[2]);
			m = Integer.parseInt(args[3]);
		} else {
        	throw new IllegalArgumentException("Incorrect number of arguments passed.");
        }
		
		
		double [] reducedTemps = new double[] { 0.6, 0.8, 1.0, 1.05, 1.1, 1.15, 1.2, 1.3, 1.4, 1.6, 2, 2.5, 3, 5, 10, 15, 20, 30, 50, 100, 500};
		//double [] reducedTemps = new double[] { 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 2, 2.5, 3, 5, 10, 15, 20, 50, 100};

		//double [] reducedTemps = new double[] { 1.0};
		
		for (int i=0;i<reducedTemps.length;i++) {
			reducedTemp = reducedTemps[i];
			
			// Number of grid points in the discretizations of r- and k-space
			
			
			// Get Mayer function for this discretization
			double[] fr = getfr( N, del_r,reducedTemp);
		
			
	        
	        
			PercusYevick py = new PercusYevick(); 
			double[] B = py.computeB(fr, m, N, del_r, false);
			
			
			
			/*
			System.out.println("Values computed here at T* = " + reducedTemp + ":");
			System.out.println("r_max = " + r_max + ", log2N = " + power);
			for (int j=2;j<=m;j++) {
				System.out.println("B" + j + " = " + B[j-2]);
			}
			*/
			
			System.out.println(reducedTemp + "  " + B[4-2]);
		}
		
		

	}

	public static double[] getfr(int N, double del_r, double reducedTemp) {
	    
		double sigma = 1.0;
		
	    double epsilon = 1.0;
		
		Space space = Space3D.getInstance();
		
		P2LennardJones p2 = new P2LennardJones(space, sigma, epsilon);
		
		
		
		double[] fr = new double[N];  // Holds discretization of Mayer function in r-space
		
		for (int n = 0; n<N; n++) {
			
			double r = n*del_r;
			
			double u = p2.u(r*r);
			
			double x = -u/reducedTemp;
			
			if ( Math.abs(x) < 0.01) {
				
				fr[n] = x + x*x/2.0 + x*x*x/6.0 + x*x*x*x/24.0 + x*x*x*x*x/120.0;	
				
			} else {
				
				fr[n] = Math.exp(x)-1.0;
				
			}
	
		}
		
		return fr;
		
	}


}
