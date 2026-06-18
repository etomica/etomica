/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial.integralequation;

import etomica.math.numerical.SineTransform;


/**
 * 
 * Calculates the Percus-Yevick (PY) virial coefficients of second to Mth order for any spherically-symmetric Mayer function, fr.
 * 
 * These values differ above third order for the compressilibility and virial routes.
 * 
 * This class has only been tested for the hard sphere and Lennard-Jones potentials.
 * 
 * @author kate
 *
 */

public class PercusYevick {

	public PercusYevick() {
	}
	
	public void setRoute(boolean compressibility) {
		this.compressibility=compressibility;
	}
	
	public double[] computeB(double[] fr, int M, int N, double del_r, boolean DCF) {
		
		/*******************************************
		/*******************************************
		 * 
		 * Computes Percus-Yevick approximations of the virial coefficients up to Mth order.  
		 * The second and third coefficients are fully accurate (the PY approximation is exact).
		 * 
		/*******************************************
		/********************************************/
	
		double[] B = new double[M];
		
		SineTransform dst = new SineTransform();
		
		double[] fk = dst.forward(fr, del_r);
		
		// Arrays to store the density expansion coefficients of c(r) and h(r)
		double[] cnr = new double[N];
		double[] hnr = new double[N];
        double[][] cnk = new double[M][0];
        double[][] hnk = new double[M][0];
		
		// Fill zeroth-order density expansion coefficients of c(r) and h(r)
		for (int i=0;i<N;i++) {
			cnr[i] = fr[i];
			hnr[i] = fr[i];
		}
		
		double B2;
		
		if (compressibility) { //virial route
			B2 = -1.0/(2.0)*(fk[0]);
		} else { //virial route
			B2 = 0;
			for (int i=0;i<N-1;i++) {
				double r = del_r*i;
				double g0 = fr[i]+1;
				B2 = B2+ g0*(rdfdr[i])*r*r;
			}
			
			double r = del_r*(N-1);
			double g0 = fr[N-1]+1;
			B2 = B2 + 0.5*g0*(rdfdr[N-1])*r*r;
			B2 = B2*4.0*Math.PI/6.0*del_r;
		}
		
		B[0] = B2;
		
		// System.out.println("B2 = " + (B2));
		
		//Compute B3 (from c1) up to BM (from c(M-2))
        cnk[0] = dst.forward(cnr, del_r);
		for (int m = 1; m <= M-2; m++) {
			
			/**************************************************************************************
			/**************************************************************************************
			 * Apply the Ornstein-Zernike relation to compute mth-order density expansion of t.
			/**************************************************************************************
			/***************************************************************************************/
			
			// we computed cnk during the last iteration
            hnk[m-1] = SineTransform.forward(hnr, del_r);
			double[] tmr = OrnsteinZernike.tCompute(cnk, hnk, N, m, del_r);
			
			/**************************************************************************************
			/**************************************************************************************
			 * Apply the Percus-Yevick approximation to compute mth-order density expansion of c.
			/**************************************************************************************
			/***************************************************************************************/
			
			for (int i = 0; i < N; i++) {
				
				cnr[i] = fr[i]*tmr[i];
				
 			}
			
			/*******************************************
			/*******************************************
			 * Update h
			/*******************************************
			/********************************************/
            
            for (int i=0; i<N; i++){
				hnr[i] = cnr[i]+ tmr[i]; 
			}	
			
			/*******************************************
			/*******************************************
			 * Calculate (m+2)th virial coefficient
			/*******************************************
			/********************************************/
			// we'll use this next iteration and here for compressibility route
			cnk[m] = dst.forward(cnr, del_r);
			
			double Bm = 0;
			if (compressibility) { //compressibility route
				Bm = -1.0/(m+2.0)*(cnk[m][0]); // B3 for m = 1
			} else { //virial route
				
				for (int i=0;i<N-1;i++) {
					double r = del_r*i;
					Bm = Bm + hnr[i]*(rdfdr[i])*r*r;
				}
				
				double r = del_r*(N-1);
				Bm = Bm + 0.5*hnr[N-1]*(rdfdr[N-1])*r*r;
				Bm = Bm*4.0*Math.PI/6.0*del_r;
			}
			
			B[m]=Bm;
            //System.out.println("B"+(m+2) + " = "+ (Bm) );
			
            
		
		}
		
		if (DCF) {
			return cnr;
		}
		return B;
	}
	
	public void setrdfdr(double[] rdfdr) {
		this.rdfdr = rdfdr; 
	}
	
	public boolean compressibility = true;
	public double[] rdfdr;
}
