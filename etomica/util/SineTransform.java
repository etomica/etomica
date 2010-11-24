package etomica.util;

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
	
	public SineTransform() {
	}
	
	
	public double[] forward(double[] f, double del_r, double r_max) {
	
		
		double del_k = Math.PI/r_max;
		//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		// Make auxiliary vector, Fr = r*fr, not including the zeroth mode (r=0)
		//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		
		int N = f.length;
		
		double[] Fr = new double[N-1]; 
		
		
		for (int n=0; n < N-1; n++) {
			
			Fr[n] = (n+1)*del_r*f[n+1];
			
		}
		
		//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		// Make auxiliary vector, Fr2 = [0,Fr,*0*,reverse(-Fr)], to evaluate DST with DFT
		//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		
		double[] Fr2 = new double[2*(N)];
		
		Fr2[0] = 0; 
		Fr2[N] = 0;
		for (int i = 0; i<N-1; i++) {
			Fr2[i+1] = Fr[i];           // from elements 1 to N-1
		    Fr2[2*N-1-i] = -Fr[i];      // from elements 2*N-1 to N+1   
		}
	
		//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		// Perform DFT on auxiliary vector Fr2
		//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

		FastFourierTransform fourier = new FastFourierTransform();
		FastFourierTransform.BACKUP=false;
		
		double[] Fr2i = new double[2*(N)]; // no imaginary part
		fourier.setData(Fr2, Fr2i);
		fourier.transform();
		
		double[] Fk2 = new double[2*(N)];
		Fk2 = fourier.getImaginary();
	
		//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		//% Extract the DST, Fk, from Fk2 = -1/N*[#,Fk,#,#]
		//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
		double[] Fk = new double[N];
		
		for (int i = 0; i<N-1; i++) {
			
			 Fk[i] = -(N)*Fk2[i+1]; 
			
		}
		
		//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		//% Compute fk (including zeroth mode) from Fk: 
		//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
	
		double [] fk = new double [N];
		
	    for (int i = 1; i<N; i++) {
		    
		    // Special consideration of zeroth mode, fk(k=0):
	    	// Even if there were an Fk(k=0) we could not utilize it in the regular way
		     
		    // fk[0] = fk[0] + 4.0*Math.PI*(Fr[i-1]*(i)*del_r)*del_r;
	    	// fk[0] = fk[0] + (4.0*Math.PI*del_r)*(Fr[i-1]*(i)*del_r);
	    	
	    	fk[0] = fk[0] + (Fr[i-1]*(i)*del_r);
	    	
	    	
		    // Modes 1 through N-1: 
	    	
		     fk[i] = 4.0*Math.PI*(Fk[i-1]/((i)*del_k))*del_r; 
	    	// fk[i] = 4.0*Math.PI*(Fk[i-1]/(i))/del_k*del_r; 
	    	// fk[i] = 4.0*Math.PI*(Fk[i-1]/(i))*(r_max/pi)*del_r; 
		    // fk[i] = 4.0*(Fk[i-1]/(i))*r_max*del_r; 
	    	
		    // Fk[0] corresponds to k=1.
		    
	    	// fk[i] = 4.0*(Fk[i-1]/(i))*r_max*del_r; // 1/del_k = r_max/pi

		}
	    
	    fk[0] = fk[0]*4.0*Math.PI*del_r;
		
		return fk;
	}
	
	public double[] reverse(double[] fk, double del_r, double r_max) {
		
		double del_k = Math.PI/r_max;
		
		//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		//% Make auxiliary vector, Fk = k*fk, ignoring the k=0 mode  
		//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		
		int N = fk.length;

	    double[] Fk = new double[N-1];
		
		for (int i=0;i<N-1;i++) {
		    
		    Fk[i] = fk[i+1]*(i+1)*del_k; 
			
			//  Fk[i] = fk[i+1]*(i+1)*Math.PI/r_max;
		    
		}
		
		//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		// Make auxiliary vector, Fk2 = [0,Fk,0,reverse(-Fk)], to evaluate DST with DFT
		//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		
		double[] Fk2 = new double[2*(N)];
	
		Fk2[0] = 0;
		Fk2[N] = 0; 
		
		for (int i=0; i<N-1; i++) {
			Fk2[i+1] = Fk[i]; 
			Fk2[2*N-1-i] = -Fk[i];    
		}
		
		//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		//% Perform DFT on auxiliary vector Fk2
		//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		
		double[] Fk2i = new double[2*(N)]; // no imaginary part
		
		FastFourierTransform fourier = new FastFourierTransform();
		FastFourierTransform.BACKUP=false;
		fourier.setData(Fk2, Fk2i);
		fourier.transform();
		
		double[] Fr2 = new double[2*(N)];
		Fr2 = fourier.getImaginary();

		//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		//% Extract Fr from Fr2
		//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		
		double[] Fr = new double[N-1];
		
		for (int i=0; i<N-1; i++) {

			  Fr[i] = -(N)*Fr2[i+1]; 

		}

		//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		//% Compute fr from Fr=r*fr;
		//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

		double[] fr = new double[N];
			
		for (int i=1; i<N; i++) {

		    //Special consideration for zeroth mode
		    fr[0] = fr[0] + 1/(2*Math.PI*Math.PI)*(Fk[i-1]*i*del_k)*del_k;
			
			//fr[0] = fr[0] + (Fk[i-1]*i); //del_k = pi/r_max

		    fr[i] = 1.0/(2.0*Math.PI*Math.PI)*(Fr[i-1]/(i*del_r))*del_k;
			
			// fr[i] = 1.0/(2.0*Math.PI)*(Fr[i-1]/(i*del_r))/r_max;  //del_k = pi/r_max
			
			//fr[i] = 1.0/(2.0*Math.PI)*(Fr[i-1]/(i*del_r))/r_max;  //del_k = pi/r_max

			
			//fr[i] = 1.0/(2.0*Math.PI)*(Fr[i]/(i))/r_max/r_max;  // 1/(r_max*del_r) = 1/(r_max*r_max/(N-1)) = (N-1)/(r_max^2)
			
			//fr[i] = 1.0/(2.0*Math.PI)*(Fr[i]/(i))/r_max/r_max;
		    
		}
		
		fr[0] = fr[0]/(2.0*r_max*r_max);
		
		return fr;
	}

	

}
