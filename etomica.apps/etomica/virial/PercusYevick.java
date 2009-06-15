package etomica.virial;

import etomica.potential.P2LennardJones;
import etomica.space.Space;
import etomica.space3d.Space3D;
import etomica.util.SineTransform;
import etomica.util.SlowDiscreteSineTransform_del_k;


/**
 * 
 * Calculates the Percus-Yevick (PY) virial coefficients of second to (m+1)th order for any spherically-symmetric Mayer function, fr.
 * 
 * This class has only been tested for the hard sphere and Lennard-Jones potentials.
 * 
 * @author kate
 *
 */

public class PercusYevick {
	
	
	
	public PercusYevick() {
	}
	
	public double[] computeB(double[] fr, int M, int N, double del_r) {
		
		/*******************************************
		/*******************************************
		 * The method relies upon Fourier transforms.  The 3-D Fourier transform simplifies to the following form upon assuming spherical symmetry:
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
		 * 
		/*******************************************
		/********************************************/
	
		double[] B = new double[M];
		
		SineTransform dst = new SineTransform();
		
		double[] dummy = new double[N];
		
		
		dummy = fr;
		double del_k = Math.PI/(del_r*N);
		double[] fk = new double[N+1];
		fk = dst.forward(dummy, del_r, del_k);
		
		// arrays to store the density expansion coefficients of the transforms of c(r) and h(r)
		double[][] cmk = new double[M][N+1];
		double[][] hmk = new double[M][N+1];
		
		cmk[0] = fk; // m = 0, NOT k = 0
		hmk[0] = fk; 
		
		double B2 = -1.0/(2.0)*(cmk[0][0]);
		
		B[0] = B2;
		
		double Bm; // virial coefficients for clusters of (m+2) molecules
		double md; // double version of m
		
		// System.out.println("B2 = " + (B2));
		
		for (int m = 1; m < M; m++) {
			
			
			/*******************************************
			/*******************************************
			 * Calculation of the mth order expansion 
			 * coefficient of the transform of t
			/*******************************************
			/********************************************/
			
			double[] tk_kvec = new double[N+1];
			double[] tr_rvec = new double[N];
			double[] cr_rvec = new double[N];
			
		
			for (int k = 0; k < N+1; k++) {
				
				double[] hk_mvec = new double[m];
				double[] ck_mvec = new double[m];
				
				for (int j=0; j<m; j++){
					
					hk_mvec[j] = hmk[j][k];
					ck_mvec[j] = cmk[j][k];
				}
				
				/*******************************************
				/*******************************************
				 * Apply the Ornstein-Zernike relation
				/*******************************************
				/********************************************/
				
				tk_kvec[k] = 0;
				
				for (int i=0; i <= m-1; i++) {
					
					int j = m-1-i;
					
					tk_kvec[k] += hk_mvec[i]*ck_mvec[j];
					
				}
					
			}
		
			
			/*******************************************
			/*******************************************
			 * Apply the Percus-Yevick approximation
			 * This must be done in r-space
			/*******************************************
			/********************************************/
			
			dummy = tk_kvec;
			tr_rvec = dst.reverse(dummy, del_r, del_k);
		    
			for (int n = 0; n < N; n++) {
				cr_rvec[n] = fr[n]*tr_rvec[n];
			}
			
			
			/*******************************************
			/*******************************************
			 * Calculate (m+2)-order virial coefficient
			/*******************************************
			/********************************************/
			
			
			dummy = cr_rvec;
			cmk[m] = dst.forward(dummy, del_r, del_k);
			
			md = m;
			
			Bm = -1.0/(md+2.0)*(cmk[m][0]); // B3 for m = 1
			
			B[m] = Bm;

            // System.out.println("B"+(m+2) + " = "+ (Bm) );
			
            /*******************************************
			/*******************************************
			 * Update h
			/*******************************************
			/********************************************/
			
			
			for (int k=0; k<N+1; k++){
				hmk[m][k] = cmk[m][k]+ tk_kvec[k]; 
			}	
		
		}
		
		return B;
	}
	
	
	

	
}
