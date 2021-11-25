/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial.integralequation;

import etomica.potential.P2HePCKLJS;
import etomica.potential.P2LennardJones;
import etomica.potential.Potential2Soft;
import etomica.space.Space;
import etomica.space3d.Space3D;
import etomica.units.Kelvin;
import etomica.util.Constants;

/**
 * 
 * Computes the Percus-Yevick compressibility-route approximation for to the 
 * first quantum correction of Bn, where the first quantum correction is 
 * as defined Kim and Henderson (1968,1966).  The gradient of u is equal
 * to d2u/dr2 + 2/r*(du/dr).
 * 
 * Considers the same potential as Kim and Henderson: Lennard-Jones parameterized
 * for helium by Haberlandt.  
 * 
 * Percus-Yevick approximation is fully accurate at second and third order, and values
 * match those of Kim and Henderson to three or four sig figs at second and third orders.
 * 
 * Applicable only for spherically symmetric pair potentials.
 * 
 * @author kate
 *
 */


public class BnPY_1QC_Helium {
	
public static void main(String[] args) {
	
		Space space = Space3D.getInstance();
		double r_max = 100; // Defines range of separation distance, r = [0 rmax]
        boolean LJ = false;
        double[] temps; double[] kTs;
        Potential2Soft p2;

        if (LJ) {
			sigma = Math.pow((10.24*1e24/Constants.AVOGADRO),1.0/3.0); //Angstroms
			System.out.println("sigma = " + sigma);
			epsilon = Kelvin.UNIT.toSim(9.66);

			temps =  new double[] {0.6, 0.7, 0.8, 0.9, 1.0, 4.0, 10, 50, 100 };
			//temps =  new double[] {0.8, 1.0, 1.2, 5, 10}; 
			//temps =  new double[] { 0.6, 0.8, 1.0, 1.2, 1.5, 2.0, 2.5, 3.0, 5.0, 10};
			kTs = new double[temps.length];
			for (int i=0;i<kTs.length;i++) {
				kTs[i] = temps[i]*epsilon;
			}
			p2 = new P2LennardJones(sigma,epsilon);
			
			//b2HS, to remove sigma dependence
			convert = (3.0/(2.0*Math.PI*sigma*sigma*sigma)); 
			
			// to include mass dependence, if desired
			//convert2 = 2.74*2.74; 
			
			//double lambda = Constants.PLANCK_H/(sigma*Math.pow(mass*epsilon, 0.5));
			//mass = Constants.PLANCK_H*Constants.PLANCK_H/(sigma*sigma*epsilon*2.74*2.74);
			//System.out.println("mass = " + mass);

        } else {
        	
        	p2 = new P2HePCKLJS();
        	temps = new double[] {83.15, 198.15, 250, 350, 400, 450, 500};
        	//temps = new double[] {350};
        	kTs = new double[temps.length]; 
			for (int i=0;i<temps.length;i++) {
				kTs[i] = Kelvin.UNIT.toSim(temps[i]);
			}
			//r_max=200;
			convert = Constants.AVOGADRO*1e-24;
			mass = 4.002602; 
			
			// to include mass dependence, if desired
			convert2 = Constants.PLANCK_H*Constants.PLANCK_H/(mass);

			
        }

		int power = 8; // Defines discretization
		int N = 1<<power;
	
		int m = 4; // highest order of virial coefficient to be calculated
		double temp;

		for (int i=0;i< temps.length;i++) {
			temp = temps[i];
			double kT = kTs[i];

			boolean printapalooza = false;
			
			if (convergence) {
				double[][] results = getConvergence(p2, power, m, r_max, kT, printapalooza);
			
				if (!printapalooza) {
					System.out.print(temp);
					
					for (int j = 0; j<m-1; j++) {	
						System.out.print("    "  + results[0][j] ); //+ "    " + results[1][j]);
					} 
					System.out.print("    2^"+results[2][0] + " grid points");
					System.out.println();
	
				}
			} else {
				
				double del_r = r_max/((double)(N-1));
				double[] B = computeB(m, N,del_r,kT,p2);
				
				System.out.print(temp + "    "  + B[1]);
				

				System.out.println();
				
			}

		}

	}

	public static double[][] getConvergence (Potential2Soft p2, int power, int m, double r_max, double kT, boolean printapalooza) {
		
		double [][] results = new double[3][m];
		double [] newBPY = new double[m];
		double [] errorBPY = new double[m];
		double [] oldBPY = new double[m];
	
		int N = (int) Math.pow(2, power);
		
		
		double error = 1.0;
		double tol = 0.0000001;
		while (error > tol) {  // increase N to improve accuracy of B2 and B3
			
			double del_r = r_max/((double)(N-1));
			
			newBPY = computeB(m, N,del_r,kT,p2);
			error=0;
			for (int i = 0; i<m-1; i++) {
				errorBPY[i] = Math.abs(newBPY[i]-oldBPY[i]);
				if (errorBPY[i]>error) {
					error = errorBPY[i];
				}
			}
	
			
			if (printapalooza) {
				System.out.print(kT);
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
	public static double[] computeB(int m, int N, double del_r, double kT, Potential2Soft p2) {

        // Get Mayer function for this discretization
        double[] fr = getfr(N, del_r, kT, p2);
        double[] rdudr = getrdudr(N, del_r, p2);
        double[] r2d2udr2 = getr2d2udr2(N, del_r, p2);

        PercusYevick1QC py = new PercusYevick1QC();
        py.setrdudr(rdudr);
        py.setr2d2udr2(r2d2udr2);

        double[] B = py.computeB(fr, m, N, del_r, false);

        double c = 1.0;
        for (int i = 0; i < m - 1; i++) {

            int n = i + 2;
            double prefactor = ((double) (n - 1)) / (24 * Math.PI * (kT / epsilon) * (kT / epsilon));
			
			c = c*convert;

			prefactor = prefactor*c*(sigma*sigma/epsilon)*convert2;

			B[i] = B[i]*prefactor;
		}
		
		return B;
		
	}

	public static double[] getfr(int N, double del_r, double kT, Potential2Soft p2) {
	    
		double[] fr = new double[N];  // Holds discretization of Mayer function in r-space
		
		for (int n = 0; n<N; n++) {
			
			double r = n*del_r;

			double u = p2.u(r*r);
			
			double x = -u/kT;
			
			if ( Math.abs(x) < 0.01) {
				
				fr[n] = x + x*x/2.0 + x*x*x/6.0 + x*x*x*x/24.0 + x*x*x*x*x/120.0;

			} else {
	
				fr[n] = Math.exp(x)-1.0;

			}
			
		}
		
		return fr;
		
	}
	
	public static double[] getrdudr(int N, double del_r, Potential2Soft p2) {
		
		double[] rdudr = new double[N];  
		
		for (int n = 1; n<N; n++) {
			
			double r = n*del_r;
			
			rdudr[n] = p2.du(r*r);

		}
		
		return rdudr;
		
	}
	
	public static double[] getr2d2udr2(int N, double del_r, Potential2Soft p2) {
		
		double[] r2d2udr2 = new double[N];  
		
		for (int n = 1; n<N; n++) {
			
			double r = n*del_r;
			
			r2d2udr2[n] = p2.d2u(r*r);

		}
		
		return r2d2udr2;
		
	}

	private static double convert = 1.0;
	private static double convert2 = 1.0;
	private static boolean convergence = true;
	private static double mass = 1.0;
	private static double sigma = 1.0;
	private static double epsilon = 1.0;

}
