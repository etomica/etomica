/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial;

import etomica.math.numerical.SineTransform;

/**
 * 
 * Computes the Percus-Yevick compressibility-route approximation to the 
 * first quantum correction of Bn, where the first quantum correction is 
 * as defined Kim and Henderson (1968,1966).  The gradient of u is equal
 * to d2u/dr2 + 2/r*(du/dr).
 * 
 * Percus-Yevick approximation is fully accurate at second and third order.
 * 
 * Applicable only for spherically symmetric pair potentials.
 * 
 * @author kate
 *
 */

public class PercusYevick1QC {

	public PercusYevick1QC() {
	}
	
	public double[] computeB(double[] fr, int M, int N, double del_r, boolean DCF) {
		
		double[] B = new double[M];
		
		SineTransform dst = new SineTransform();
		double r_max = del_r*(N-1);
		
		double[] dummy = new double[N];
		dummy = fr;
		double[] fk = new double[N];
		fk = dst.forward(dummy, del_r);
		
		// Arrays to store the density expansion coefficients of c(r) and h(r)
		double[][] cnr = new double[M][N];
		double[][] hnr = new double[M][N];
		
		// Fill zeroth-order density expansion coefficients of c(r) and h(r)
		for (int i=0;i<N;i++) {
			cnr[0][i] = fr[i];
			hnr[0][i] = fr[i];
		}
		
		double B2 = 0;


		for (int i=0;i<N-1;i++) {
			double e12 = fr[i]+1.0;
			B2 = B2 + e12*(r2d2udr2[i] + 2.0*rdudr[i])*del_r;
		}
		double e12 = fr[N-1]+1.0;
		B2 = B2 + 0.5*e12*(r2d2udr2[N-1] + 2.0*rdudr[N-1])*del_r;

		
		B[0] = B2;
		
		// System.out.println("B2 = " + (B2));
		double[] cmr = new double[N];
		
		//Compute B3 (from c1) up to BM (from c(M-2))
		for (int m = 1; m <= M-2; m++) {
			
			/**************************************************************************************
			/**************************************************************************************
			 * Apply the Ornstein-Zernike relation to compute mth-order density expansion of t.
			/**************************************************************************************
			/***************************************************************************************/
			
			OrnsteinZernike oz = new OrnsteinZernike();
			
			double[] tmr = oz.tCompute(cnr, hnr, N, m, del_r);
			
			/**************************************************************************************
			/**************************************************************************************
			 * Apply the Percus-Yevick approximation to compute mth-order density expansion of c.
			/**************************************************************************************
			/***************************************************************************************/
			
			for (int i = 0; i < N; i++) {
				
				cnr[m][i] = fr[i]*tmr[i];
				
				cmr[i] = cnr[m][i];
				
 			}
			
			/*******************************************
			/*******************************************
			 * Update h
			/*******************************************
			/********************************************/
            
            for (int i=0; i<N; i++){
				hnr[m][i] = cnr[m][i]+ tmr[i]; 
			}	
			
			/*******************************************
			/*******************************************
			 * Calculate (m+2)th virial coefficient
			/*******************************************
			/********************************************/
			
			dummy = cmr;
			double[] cmk = new double[N];
			cmk = dst.forward(dummy, del_r);
			
			double Bm = 0;

			for (int i=0;i<N-1;i++) {
				Bm = Bm + (hnr[m][i])*(r2d2udr2[i]+2.0*rdudr[i])*del_r;
			}
			Bm = Bm + 0.5*(hnr[m][N-1])*(r2d2udr2[N-1]+2.0*rdudr[N-1])*del_r;
			
			B[m] = Bm;

		}
		
		return B;

	}

	public void setrdudr(double[] rdudr) {
		this.rdudr = rdudr; 
	}
	
	public void setr2d2udr2(double[] r2d2udr2) {
		this.r2d2udr2 = r2d2udr2; 
	}
	
	public double[] rdudr;
	public double[] r2d2udr2;
}

