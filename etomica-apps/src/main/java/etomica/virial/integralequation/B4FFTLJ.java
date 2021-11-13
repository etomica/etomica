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
 * Computes four approximations of B4 by FFT: HNC(C), HNC(V), PY(C), and PY(V).
 * 
 * HNC: hypernetted-chain integral-equation theory
 * PY: Percus-Yevick integral-equation theory
 * (C):compressibility route 
 * (V):virial (pressure) route 
 *  
 * The HNC(V) approximation is equivalent to the "simple part" of the f-bond-only formulation for (at least) B4, B5, and B6.
 * 
 * Values agree with those computed through the recursive formalism (the Ornstein-Zernike formulation) in BnHNCLJ.java and BnPYLJ.java.
 * 
 * Diagrams are labeled using the nomenclature of Barker, Leonard, and Pompe (1966): 
 * "Fifth Virial Coefficients", JCP 44(11); e.g., D4 is the four-membered ring.
 * 
 * X is added to the diagram name to denote if an x=r*df/dr bond is present.  This our own nomenclature...
 * 
 * @author Kate Shaul
 *
 */


public class B4FFTLJ {
	
public static void main(String[] args) {
	
	double [] temps = new double[] { 0.6, 0.8, 1.0, 1.05, 1.1, 1.15, 1.2, 1.3, 1.4, 1.6, 2, 2.5, 3, 5, 10, 15, 20, 30, 50, 100, 500};
	
	for (int t=0; t<temps.length; t++) {
		
		double reducedTemp = temps[t]; // kT/epsilon
		
		int power = 13; // Defines discretization
		double r_max = 100; // Defines range of separation distance, r = [0 rmax]
		int N = 1<<power;  // Number of grid points in the discretizations of r- and k-space

		double del_r = r_max/(N-1);
		
		double[] fr = getfr( N, del_r,reducedTemp);

		
		SineTransform dst = new SineTransform();
		
		double[] fk = dst.forward(fr, del_r);
		
		double[] ffk = new double[N];
		double[] fffk = new double[N];
		for (int i = 0;i<N; i++) {
			ffk[i] = fk[i]*fk[i];
			fffk[i] = fk[i]*fk[i]*fk[i]; 
		}
		
		double[] ffr = dst.reverse(ffk, del_r);
		double[] d1r = dst.reverse(fffk, del_r);

		
		Space space = Space3D.getInstance();	
		P2LennardJones p2 = new P2LennardJones(space, 1, 1);
		
		double C3 = 0;
		double D4 = 0;
		double D5 = 0;
		double D5X = 0;
		
		for (int i = 1;i<(N-1); i++) {
			
			double c = 4.0*Math.PI*(i*del_r)*(i*del_r)*del_r;
			
			C3 = C3 + fr[i]*ffr[i]*c;
			 
			D4 = D4 + fr[i]*d1r[i]*c;
			
			D5 = D5 + fr[i]*(ffr[i]*ffr[i])*c;
			
			double r = del_r*i;
			double x = (fr[i]+1)*(-1.0/reducedTemp)*p2.du(r*r);
			D5X = D5X + x*(ffr[i]*ffr[i])*c;
			
		}
		
		double c = 4.0*Math.PI*r_max*r_max*del_r;
		
	    C3 = C3 + 0.5*fr[N-1]*ffr[N-1]*c;
		
		D4 = D4 + 0.5*fr[N-1]*d1r[N-1]*c; 
		
		D5 = D5 + 0.5*fr[N-1]*( ffr[N-1]*ffr[N-1] )*c; 
		
		double r = del_r*(N-1);
		double x = (fr[N-1]+1)*(-1.0/reducedTemp)*p2.du(r*r);
		D5X = D5X + x*(ffr[N-1]*ffr[N-1])*c;
		
		double B3 = -1.0/(3.0)*C3;
		
		double B4HNCV = -3.0/8.0*( D4 + 2.0*D5);
		
		double B4HNCC = -3.0/8.0*D4 + -5.0/8.0*D5;
		
		double B4PYV  = -3.0/8.0*( D4 + 2.0*D5)-1.0/12.0*D5X;
			
		double B4PYC  = -1.0/4.0*( D4 + 2.0*D5);
		
		System.out.println(reducedTemp + "   "  + (B4HNCC)+ "   "  + (B4HNCV)+ "   "  + (B4PYC)+ "   "  + (B4PYV));
		
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
