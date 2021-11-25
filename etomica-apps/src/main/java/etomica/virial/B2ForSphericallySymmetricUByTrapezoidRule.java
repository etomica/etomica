/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial;

import etomica.potential.P2ArgonAziz1993;
import etomica.potential.Potential2Soft;
import etomica.space.Space;
import etomica.space3d.Space3D;

/* 
 * As the name implies, this class computes the second virial coefficient for a spherically symmetric potential using the trapezoid rule.  
 * 
 * Kate Shaul, 2010
 */

public class B2ForSphericallySymmetricUByTrapezoidRule {
	
	public B2ForSphericallySymmetricUByTrapezoidRule() {
		
	}
	
	public static void main(String[] args) {
		
		double r_max = 100.0;
		
		System.out.println("Argon pair potential of Aziz 1993");
		System.out.println();
		System.out.println("T (K)    B2       difference between value and that computed with half as many points");
		
		Space space = Space3D.getInstance();
		P2ArgonAziz1993 p2 = new P2ArgonAziz1993();
		
		int power = 10;
		
		double[] temps = { 100, 133.15, 150, 200, 250, 300, 400, 500, 600, 700, 800, 900, 1000 }; // Kelvin
		
		
		boolean printapalooza = false;
		
		for (int t=0; t < temps.length; t++) {
			
		
			double temp = temps[t]; // kT/epsilon
			
			boolean allFFT = false;
			
			if (args.length == 0) {
				
			} else if (args.length == 3) { 
				
				temp = Double.parseDouble(args[0]);
				power = Integer.parseInt(args[1]);
				r_max = Double.parseDouble(args[2]);
				
				
			} else if (args.length == 6) {
				temp = Double.parseDouble(args[0]);
				power = Integer.parseInt(args[1]);
				r_max = Double.parseDouble(args[2]);
				
				allFFT = Boolean.parseBoolean(args[5]);
				
			} else {
	        	throw new IllegalArgumentException("Incorrect number of arguments passed");
	        }
		
	
		
		    double[] results = getResults (p2, power, r_max, temp, printapalooza);
		
		    double B2 = results[0];
		    double error = results[1];
		    double power2 = results[2];
			
			System.out.println(temp + "     " + B2 + "    " + error + "    " + power2);
			
		}


	}

    public static double[] getResults (Potential2Soft p2, double power, double r_max, double temp, boolean printapalooza) {
    	
    	double B2 = 0;
		int count = 0;
		double[] powers = new double[30];
		double[] errors = new double[30];
		double power2 = power;
		double error = 1.0;
		double [] results = new double[3];
		while (error > 1e-3) {
		
			// Number of grid points in the discretization of r-space;
	       
			int N = (int) Math.pow(2, power2) - 1;
			double del_r = r_max/(N-1);
			
			double[] fr = getfr( p2, N, del_r,temp);
			 
			double newB2 = 0;
			
			// Trapezoid rule 
			for (int i=1;i<N-1; i++) {	 
				 
				double r = i*del_r;
				 
				 newB2 += r*r*fr[i]; 
				
			}
			
			newB2 += 0.5*del_r*del_r*(N-1)*(N-1)*fr[N-1];
			
			newB2 = -2.0*Math.PI*del_r*newB2*0.60221415; 
			
			error = Math.abs(B2-newB2);
			
			powers[count] = power2;
			errors[count] = error;

		
			B2 = newB2 ;
			
			if (printapalooza) {
			
				System.out.println(power2 + "   " + B2 + "  " + error);
			
			}
			
			results[0] = B2; 
			results[1] = error;
			results[2] = power2;
			
			count++;
			power2 = power2+1.0;
			
		}
		
		
		return results;
		
    }

	public static double[] getfr(Potential2Soft p2, int N, double del_r, double temp) {
		
		double[] fr = new double[N];  // Holds discretization of Mayer function in r-space
		
		for (int n = 0; n<N; n++) {
			
			double r = n*del_r;
			
			
			
			//System.out.println(r + "  " + x);
			
			if ( r < 2.5) {
				
				fr[n]=-1.0;
				
			} else {
				
				double u = p2.u(r*r); // returned in Kelvin

				double x = -u/temp;
				
				if ( Math.abs(x) < 0.01) {
					
					fr[n] = x + x*x/2.0 + x*x*x/6.0 + x*x*x*x/24.0 + x*x*x*x*x/120.0;
					
					
				} else {
					
					
					
					fr[n] = Math.exp(x)-1.0;
					
				}
				
			}
			
		}
		
		return fr;
		
	}

}
