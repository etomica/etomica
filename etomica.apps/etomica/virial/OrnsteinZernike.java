package etomica.virial;

import etomica.util.SineTransform;


/**
 * @author kate
 *
 */

public class OrnsteinZernike {
	
	
	
	public OrnsteinZernike() {
	}
	
	public double[] tCompute(double[][] cnr, double[][] hnr, int N, int m, double del_r) {
		
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

		// Create transforms hnk and cnk
		
		double[][] hnk = new double[m][N];
		double[][] cnk = new double[m][N];
	
		double[] cnrCopy = new double[N];
		double[] hnrCopy = new double[N];
		double[] cnkCopy = new double[N];
		double[] hnkCopy = new double[N];
		
		SineTransform dst = new SineTransform();
		
		double r_max = del_r*(N-1);
		
		for (int n=0; n<(m); n++) {
			
			for (int i=0; i<N; i++) {
			
				cnrCopy[i] = cnr[n][i];
				
				hnrCopy[i] = hnr[n][i];
			}
			
			cnkCopy = dst.forward(cnrCopy, del_r);
			
			hnkCopy = dst.forward(hnrCopy, del_r);
		
			for (int i=0; i<N; i++) {
				
				cnk[n][i] = cnkCopy[i];
				
				hnk[n][i] = hnkCopy[i];
				
			}
		
		}

		double[] tmk = new double[N];
		
		for (int k = 0; k < N; k++) {
			
			double[] hkn = new double[m];
			double[] ckn = new double[m];
			
			for (int j=0; j<m; j++){
				
				hkn[j] = hnk[j][k];
				ckn[j] = cnk[j][k];
			}
				
			tmk[k] = 0;
			
			for (int i=0; i <= m-1; i++) {
					
				int j = m-1-i;
				
				tmk[k] += hkn[i]*ckn[j];
				
			}
				
		}
		
		double[] tm = new double[N];
		tm = dst.reverse(tmk, del_r);
	    
		return tm;	
		
	}
	
}
