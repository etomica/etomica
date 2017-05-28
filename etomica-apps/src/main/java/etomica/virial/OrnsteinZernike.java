/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial;

import etomica.math.numerical.SineTransform;


/**
 * @author kate
 *
 */

public class OrnsteinZernike {
	
	public static double[] tCompute(double[][] cnk, double[][] hnk, int N, int m, double del_r) {
		
		/*******************************************
		/*******************************************
		 * 
		 * Given lower-order density-expansion coefficients of the direct and total correlation functions (cnr and hnr arrays, respectively),  
		 * this method computes the mth-order expansion coefficient of indirect correlation function, tmr, via the Ornstein-Zernike relation. 
		 * 
		 * The Ornstein-Zernike relation is a convolution integral, best handled in Fourier space.  
		 * The 3-D Fourier transform simplifies to the following form upon assuming spherical symmetry:
		 * 
		 * f(k) = integral over (constant*f(r)*r*sin(k*r)) dr
		 * 
		 * Following Duh and Haymet (1995), the transform becomes an an even simpler sine transform if one uses auxiliary functions of the form, F(r) = r*f(r)
		 * 
		 * f(k) = integral over (constant*F(r)*sin(k*r))dr
		 * 
		 * The limit as k-->0 or r-->0 is used to evaluated the zeroth modes.
		 * 
		 * The requisite transforms are carried out by SineTransform.java and FastFourierTrasnform.java.
		 *  - The 1-D sine transform is performed using a 1-D FFT algorithm.
		 * 
		/*******************************************
		/********************************************/

		double[] tmk = new double[N];
		
		for (int k = 0; k < N; k++) {
			
			tmk[k] = 0;
		}
			
		for (int i=0; i <= m-1; i++) {
				
			int j = m-1-i;
	        for (int k = 0; k < N; k++) {
			
	            tmk[k] += hnk[i][k]*cnk[j][k];
	        }
			
		}
			
		double[] tm = SineTransform.reverse(tmk, del_r);
	    
		return tm;	
		
	}
	
}
