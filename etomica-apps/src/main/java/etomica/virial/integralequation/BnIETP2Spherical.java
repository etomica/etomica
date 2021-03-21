/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial.integralequation;

import etomica.potential.*;
import etomica.space.Space;
import etomica.space3d.Space3D;
import etomica.units.Kelvin;
import etomica.util.Constants;
import etomica.virial.PercusYevick;

/**
 * 
 * This is just a main method for HypernettedChain.java or PercusYevick.java...
 * 
 * It creates a discretization of the Mayer function for a spherically symmetric potential, which HypernettedChain (or PercusYevick) then employs to compute 
 * approximations of the virial coefficients up to mth order.  The second and third coefficients are fully accurate (the 
 * HNC (or PY) approximation is exact).
 * 
 * Use the boolean compressibilityRoute to select whether the compressibility-route or virial-route approximation is used.
 * Use the boolean QFH to select whether the quadratic Feynman-Hibbs modification to the potential is used.
 * 
 * @author kate
 *
 */


public class BnIETP2Spherical {
	
public static void main(String[] args) {
	
		Space space = Space3D.getInstance();
		double r_max = 200; // Defines range of separation distance, r = [0 rmax]
        boolean LJ = false;
        double[] temps;
        double[] simTemps;
        Potential2SoftSpherical p2;
		if (LJ) {
		// To make sure that everything is working fine:
		System.out.println("Literature values for LJ with sigma = 1 and T* = 1");
		System.out.println("B2 = -5.3158 (Sun & Teja 1996)");
		System.out.println("B3 =  1.8849 (Sun & Teja 1996) ");
		System.out.println("B4PY(c) = -2.9394 (Dyer et al 2001)");
		System.out.println("B4PY(c) = -2.925 (imprecise value from Henderson, Kim, and Oden (1966))");
		System.out.println("B4PY(v) = -6.627 (imprecise value from Henderson, Kim, and Oden (1966))");
		System.out.println("B5PY(c) = -37.402 (Barker et al 1966)\n");

		p2 = new P2LennardJones(space, 1.0, 1.0);
		//temps = new double[] { 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 2, 2.5, 3, 5, 10, 15, 20, 50, 100};
		
		temps = new double[] { 0.6, 0.8, 1.0, 1.05, 1.1, 1.15, 1.2, 1.3, 1.4, 1.6, 2, 2.5, 3, 5, 10, 15, 20, 30, 50, 100, 500};

		simTemps = temps;
		//double [] temps = new double[] { 1.0};
		
		} else {
		
			p2 = new P2HePCKLJS(space);
		
			if (PercusYevick) {
				temps = new double[] {20.0,24.5561,30.0,40.0,50.0,63.15,75.0,83.15,98.15,103.15,113.15,123.15,143.15,148.15,158.15,173.15,183.15,198.15,223.15,248.15,253.15,273.15,273.16,293.16,298.15,323.15,323.16,348.15,350,373.15,375,398.15,400,423.15,425,450,475,500};
			} else {
				temps = new double[] {20.0,24.5561,30.0,40.0,50.0,63.15,75.0};
			}
			
			//temps = new double[] {250,350,400,450};
			simTemps = new double[temps.length]; 
			for (int i=0;i<temps.length;i++) {
				simTemps[i] = Kelvin.UNIT.toSim(temps[i]);
			}
			
			mass = 4.002602;
			convert = Constants.AVOGADRO*1e-24;
		}
		int power = 15; // Defines discretization
		int N = 1<<power;

		
		
		int m = 5; // highest order of virial coefficient to be calculated
		double temp;
		if (args.length == 0) {
		}
		else if (args.length == 4) {
			power = Integer.parseInt(args[0]);
			r_max = Double.parseDouble(args[1]);
			temp = Double.parseDouble(args[2]);
			m = Integer.parseInt(args[3]);
		} else {
        	throw new IllegalArgumentException("Incorrect number of arguments passed.");
        }

		
		
		
		
		for (int i=0;i< temps.length;i++) {
			temp = temps[i];
			double simTemp = simTemps[i];
			
			// Number of grid points in the discretizations of r- and k-space
			
			
			boolean printapalooza = false;
			
			if (convergence) {
				double[][] results = getConvergence(p2, power, m, r_max, simTemp, printapalooza);
				
			
			
				if (!printapalooza) {
					System.out.print(temp);
					
					for (int j = 0; j<m-1; j++) {	
						System.out.print("    "  + results[0][j] + "    " + results[1][j]);
					}
					System.out.print("    2^"+results[2][0] + " grid points");
					System.out.println();
					
					
				}
			} else {
				
				double del_r = r_max/((double)(N-1));
				double[] B = computeB(m, N,del_r,simTemp,p2);
				
				System.out.print(temp + "    "  + B[2]);
				

				System.out.println();
				
			}

			
		}
		
		

	}

	public static double[][] getConvergence (Potential2SoftSpherical p2, int power, int m, double r_max, double temp, boolean printapalooza) {
	
		double [][] results = new double[3][m];
		double [] newBPY = new double[m];
		double [] errorBPY = new double[m];
		double [] oldBPY = new double[m];
	
		int N = (int) Math.pow(2, power);
		
		
		double error = 1.0;
		double tol = 0.0001;
		while (error > tol) {  // increase N to improve accuracy of B2 and B3
			
			double del_r = r_max/((double)(N-1));
			
			newBPY = computeB(m, N,del_r,temp,p2);

			error=0;
			for (int i = 0; i<m-1; i++) {
				
				newBPY[i] = newBPY[i];
				
				errorBPY[i] = Math.abs(newBPY[i]-oldBPY[i]);
				if (errorBPY[i]>error) {
					error = errorBPY[i];
				}
			}

			
			if (printapalooza) {
				System.out.print(temp);
				for (int i = 0; i<m-1; i++) {	
					System.out.print("    "  + newBPY[i] + "    " + errorBPY[i]);
				}
				System.out.println();
			}
			
				
			N=2*N;
			power=power+1;
			for (int i = 0; i<m-1; i++) {	
				oldBPY[i] = newBPY[i];
			}
				
		}
		
		for (int i = 0; i<m-1; i++) {	
			results[0][i] = newBPY[i];
			results[1][i] = errorBPY[i];
		}
		
		results[2][0] = power;
		
		return results;
	}


	public static double[] computeB(int m, int N, double del_r, double temp, Potential2SoftSpherical p2) {
		
		// Get Mayer function for this discretization
		
		double[] fr = new double[N]; double[] rdfdr = new double[N];
		if (QFH) {
			P2EffectiveFeynmanHibbs p2E = new P2EffectiveFeynmanHibbs(Space3D.getInstance(), p2);
			p2E.setTemperature(temp);		
			p2E.setMass(mass);
			fr = getfr( N, del_r,temp, p2E);
			
			
		} else {
			fr = getfr( N, del_r,temp, p2);
			rdfdr = getrdfdr( N, del_r,temp,p2);
		}

		
		double[] B;
        if (PercusYevick) {
			etomica.virial.PercusYevick py = new PercusYevick();
			py.setRoute(compressibilityRoute);

			if (!compressibilityRoute && QFH) {
				//cannot yet do virial route yet for P2EffectiveFeynmanHibbs because gradient method not present
				throw new RuntimeException("Cannot do this yet");
			}

			py.setrdfdr(rdfdr);
			B = py.computeB(fr, m, N, del_r, false);	

        } else {
			HypernettedChain hnc = new HypernettedChain();


			hnc.setRoute(compressibilityRoute);
			hnc.setrdfdr(rdfdr);

			B = hnc.computeB(fr, m, N, del_r, false);
		}
		
		
		
		double c = 1.0;
		for (int i = 0; i<m-1; i++) {
		
			c = c*convert;

			B[i] = B[i]*c;
		}
		
		return B;
		
	}

	public static double[] getfr(int N, double del_r, double temp, Potential2Spherical p2) {
	    
		double[] fr = new double[N];  // Holds discretization of Mayer function in r-space
		
		for (int n = 0; n<N; n++) {
			
			double r = n*del_r;
			
			double u = p2.u(r*r);
			
			double x = -u/temp;
			
			if ( Math.abs(x) < 0.01) {
				
				fr[n] = x + x*x/2.0 + x*x*x/6.0 + x*x*x*x/24.0 + x*x*x*x*x/120.0;	
				
			} else {
				
				fr[n] = Math.exp(x)-1.0;
				
			}
	
		}
		
		return fr;
		
	}
	
	public static double[] getrdfdr(int N, double del_r, double temp, Potential2SoftSpherical p2) {		
		
		double[] rdfdr = new double[N];  // Holds discretization of Mayer function in r-space
		
		for (int n = 1; n<N; n++) {
			
			double r = n*del_r;
			
			//To match Henderson and Oden (1966):
			//rdfdr[n] = (-Math.exp(x)/reducedTemp)*p2.du(r*r);
			// My gn and hn include the factor of e; above would double count it:
			rdfdr[n] = (-1.0/temp)*p2.du(r*r);


		}
		
		return rdfdr;
		
	}
	
	private static boolean convergence = false;
	private static boolean compressibilityRoute = false;
	private static boolean PercusYevick = false;
	private static boolean QFH = false;
	private static double convert = 1.0;
	private static double mass = 1.0;


}
