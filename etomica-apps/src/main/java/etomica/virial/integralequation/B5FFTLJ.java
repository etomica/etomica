/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial.integralequation;

import etomica.math.numerical.SineTransform;
import etomica.potential.P2LennardJones;
import etomica.space.Space;
import etomica.space3d.Space3D;

/**
 * 
 * Computes four approximations of B5 by FFT: HNC(C), HNC(V), PY(C), and PY(V).
 * 
 * HNC: hypernetted-chain integral-equation theory
 * PY: Percus-Yevick integral-equation theory
 * (C):compressibility route 
 * (V):virial (pressure) route 
 *  
 * The HNC(V) approximation is equivalent to the "simple part" of the f-bond-only formulation of B5.
 * 
 * Expressions and values for individual diagrams agree with those of Kim, Henderson, and Oden (1966):
 * "Theory of Fluids and the Fifth Virial Coefficient", JCP 45(11).
 * 
 * Note that the expressions for B5HNC(C) and B5HNC(V) are reversed in the paper.
 * 
 * Diagrams are labeled using the nomenclature of Barker, Leonard, and Pompe (1966): 
 * "Fifth Virial Coefficients", JCP 44(11); e.g., E5 is the five-membered ring.
 * 
 * X is added to the diagram name to denote if an x=r*df/dr bond is present.  This our own nomenclature...
 * 
 * @author Kate and Srihari
 *
 */


public class B5FFTLJ {
	
public static void main(String[] args) {
	
	double [] temps = new double[] { 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 2, 2.5, 3, 5, 10, 15, 20, 50, 100};
	
	for (int t=0; t<temps.length; t++) {
		
		double reducedTemp = temps[t]; // kT/epsilon
		
		int power = 17; // Defines discretization
		double r_max = 100; // Defines range of separation distance, r = [0 rmax]
		//int N = (int) Math.pow(2, power);  // Number of grid points in the discretizations of r- and k-space
		int N = 1<<power;
	
        double del_r = r_max/(N-1);
		
		double[] fr = getfr( N, del_r,reducedTemp);
		double[] dummy = new double[N];
		for (int i = 0;i<N; i++) {
			dummy[i] = fr[i];
		}
		SineTransform dst = new SineTransform();

		double[] fk = dst.forward(fr, del_r);
		
		double[] ffk = new double[N];
		double[] fffk = new double[N];
		double[] ffffk = new double[N];
		double[] fck = new double[N];
		
		for (int i = 0;i<N; i++) {
			ffk[i] = fk[i]*fk[i];
			fffk[i] = fk[i]*fk[i]*fk[i];
			ffffk[i] = fk[i]*fk[i]*fk[i]*fk[i];
		}
		
		double[] ffffr = dst.reverse(ffffk, del_r);
		double[] fffr  = dst.reverse(fffk, del_r);
		double[] ffr   = dst.reverse(ffk, del_r);
		
		double[] E6ar = new double[N];
		double[] E6br = new double[N];
		double[] E7gr = new double[N];
		double[] E7acr = new double[N];	
		
		for (int i = 1;i<(N-1); i++) {
			
			E6ar[i] = ffr[i]*fr[i]*fffr[i];
			E6br[i] = ffr[i]*ffr[i]*ffr[i];
			E7gr[i] = ffr[i]*ffr[i]*ffr[i]*fr[i];
			E7acr[i] = fr[i]*ffr[i];

		}
		
		double[] E7ack = dst.forward(E7acr, del_r);
		double[] E7ac1k = new double[N];
		
		for (int i = 0;i<N; i++) {
			E7ac1k[i] = fk[i]*E7ack[i];
		}
		
		double[] E7ac1r = dst.reverse(E7ac1k, del_r);
		double[] E7ar = new double[N];
		
		Space space = Space3D.getInstance();
		P2LennardJones p2 = new P2LennardJones(space, 1, 1);
		
		double[] E6axr = new double[N];
		double[] E7axr = new double[N];
		double[] E5xr = new double [N];
		
		for (int i = 1;i<(N-1); i++) {
			
			E7ar[i] = ffr[i]*fr[i]*E7ac1r[i];
			
			double r=i*del_r;
			double x = (fr[i]+1.0)*(-1.0/reducedTemp)*p2.du(r*r);
			
			E5xr[i] = ffffr[i]*x;
				
			E6axr[i] = ffr[i]*x*fffr[i];
				
			E7axr[i] = ffr[i]*x*E7ac1r[i];
		}
		
		// Trapezoid rule for final layer of integration
		
		double E5 = 0;
		double E6B = 0;
		double E6A = 0;
		double E7A =0;
		double E7G = 0;
		
		double E5X = 0;
		double  E6AX = 0;
		double  E7AX = 0;

		for (int i = 1;i<(N-1); i++) {
			
			double c = 4.0*Math.PI*(i*del_r)*(i*del_r)*del_r;
						 
			E5 = E5 + fffr[i]*ffr[i]*c;
			
			E6B = E6B + E6br[i]*c;
		
			E6A = E6A + E6ar[i]*c;
			
			E7A = E7A + E7ar[i]*c;
			
			E7G = E7G + E7gr[i]*c;
			
			E5X = E5X + E5xr[i]*c;
			
			E6AX = E6AX + E6axr[i]*c;
	    	
	    	E7AX = E7AX + E7axr[i]*c;

		}
		
		double c = 4.0*Math.PI*(N-1)*del_r*(N-1)*del_r*del_r;
		
    	E5 = E5 + 0.5*( fffr[N-1]*ffr[N-1] )*c; 
    	
    	E6B = E6B + 0.5*( E6br[N-1] )*c;
    	
    	E6A = E6A + 0.5*( E6ar[N-1] )*c;
    	
    	E7A = E7A + 0.5*( E7ar[N-1] )*c;
    	
    	E7G = E7G + 0.5*( E7gr[N-1] )*c;
    	
    	E5X = E5X + 0.5*(E5xr[N-1])*c;
    	
    	E6AX = E6AX + 0.5*(E6axr[N-1])*c;
    	
    	E7AX = E7AX + 0.5*(E7axr[N-1])*c;


    	double B5HNCC = -2.0/5.0*E5 + -8.0/5.0*E6A + -7.0/5.0*E7A + -7.0/30.0*(E6B+E7G);
    	
		double B5HNCV = (-2.0/5.0*E5 + -2.0*(E6A + E7A) + -1.0/(3.0)*(E6B+E7G));
		
		double B5PYC = (-1.0/5.0*E5 - E6A - E7A);
		
		double B5PYV = -2.0/5.0*E5 + -2.0*(E6A + E7A) + -1.0/(6.0)*(E6AX + 2.0*E7AX);

		System.out.println(temps[t] + "   "  + (B5HNCC)+ "   "  + (B5HNCV)+ "   "  + (B5PYC)+ "   "  + (B5PYV));
	
		//Checks:
		//double b = Standard.B2HS(1.0);
		//System.out.println(temps[t] + "   " + (E5/(b*b*b*b)) + "   " + (1.0/6.0*E5X/(-2.0/5.0)/(b*b*b*b)) + "   " + (B5PYV));
		//System.out.println(temps[t] + "   " + (E6A/(b*b*b*b)) + "   " + (E7A/(b*b*b*b)) );


	}
	
	
	}

           

	public static double[] getfr(int N, double del_r, double reducedTemp) {
		
		// Lennard-Jones Potential
	    double sigma = 1.0;
		double epsilon = 1.0;
		Space space = Space3D.getInstance();
		P2LennardJones p2 = new P2LennardJones(space, sigma, epsilon);
		double u;
		
	
		
		double[] fr = new double[N];  // Holds discretization of Mayer function in r-space
		
		for (int n = 0; n<N; n++) {
			
			double r = n*del_r;
			
			u = p2.u(r*r);
			
			double x = -u/reducedTemp;
			
			if ( Math.abs(x) < 0.01) {
				
				fr[n] = x + x*x/2.0 + x*x*x/6.0 + x*x*x*x/24.0 + x*x*x*x*x/120.0;
				
				
			} else {
				
				fr[n] = Math.exp(x)-1.0;
				
			}
			
			// fr[n] = Math.exp(x)-1.0;
			
	
		}
		
		return fr;
		
	}


}
