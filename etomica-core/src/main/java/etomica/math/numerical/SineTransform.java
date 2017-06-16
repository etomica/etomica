/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.math.numerical;

/**
 * 3D Fourier transforms of a function, f, simplify to 1D sine transforms of the auxiliary function F(r) = r*f(r) 
 * when f is spherically symmetric.  
 * 
 * Discrete sine transforms can be carried out by a fast Fourier transform operating on the vector [0 F 0 reverse(-F)].  
 * FastFourierTransform.java is used as the FFT black box.
 *  
 * @author kate
 *
 */

public class SineTransform {

	public static double[] forward(double[] f, double del_r) {
	
		
		//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		// Make auxiliary vector, Fr = r*fr
		//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		
		int N = f.length;

		//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		// Make auxiliary vector, Fr2 = [0,Fr[1:N-1],*0*,reverse(-Fr[1:N-1])], to evaluate DST with DFT
		//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		
        double[] Fr2 = new double[2*N];
		
		for (int i = 1; i<N; i++) {
			Fr2[i] = (i)*del_r*f[i];
			Fr2[N+i] = -(N-i)*del_r*f[N-i];
		}
		
		
		
		//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		// Perform DFT on auxiliary vector Fr2
		//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

		FastFourierTransform fourier = new FastFourierTransform();
		
		double[] Fr2i = new double[2*N]; // no imaginary part
		fourier.setData(Fr2, Fr2i);
		fourier.transform();
		
		double[] Fk2 = fourier.getImaginary();

		//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		//% Compute fk (including zeroth mode) from Fk: 
		//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
        double del_k = Math.PI/((N)*del_r);
	
		double [] fk = new double [N];
		
	    for (int i = 1; i<N; i++) {
		    
		    // Special consideration of zeroth mode, fk(k=0):
	    	// Even if there were an Fk(k=0) we could not utilize it in the regular way
	    	fk[0] = fk[0] + 4.0*Math.PI*(i*del_r*f[i]*(i)*del_r)*del_r;


	    	// Modes 1 through N-1: 
	    	
            fk[i] = 4.0*Math.PI*(-(N)*Fk2[i]/((i)*del_k))*del_r; 
		    // Fk[0] corresponds to k=1.

		}
	    
		
		return fk;
	}
	
	public static double[] reverse(double[] fk, double del_r) {
		
		//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		//% Make auxiliary vector, Fk = k*fk, ignoring the k=0 mode  
		//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		
		int N = fk.length;
		
		double del_k = Math.PI/((N)*del_r);
		
		//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		// Make auxiliary vector, Fk2 = [0,Fk[1:N-1],0,reverse(-Fk[1:N-1])], to evaluate DST with DFT
		//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		
		double[] Fk2 = new double[2*N];
	
		Fk2[0] = 0;
		Fk2[N] = 0; 
		
		for (int i=1; i<N; i++) {
			Fk2[i] = fk[i]*(i)*del_k;
			Fk2[2*N-i] = -Fk2[i];    
		}
		
		//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		//% Perform DFT on auxiliary vector Fk2
		//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		
		double[] Fk2i = new double[2*(N)]; // no imaginary part
		
		FastFourierTransform fourier = new FastFourierTransform();
		fourier.setData(Fk2, Fk2i);
		fourier.transform();
		
		double[] Fr2 = fourier.getImaginary(); // = Fk2i

		//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		//% Compute fr from Fr=r*fr;
		//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

		double[] fr = new double[N];
			
		for (int i=1; i<N; i++) {

		    //Special consideration for zeroth mode
            fr[0] = fr[0] + 1.0/(2.0*Math.PI*Math.PI)*(fk[i]*(i)*del_k*(i)*del_k)*del_k;

            fr[i] = 1.0/(2.0*Math.PI*Math.PI)*(-(N)*Fr2[i]/(i*del_r))*del_k;
            
            
		}
		
		
		return fr;
	}

	

}
