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
	
	
	public double[] forward(double[] f, double del_r, double del_k) {
	
		//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		// Make auxiliary vector, Fr = r*fr
		//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		
		int N = f.length;
		
		double[] Fr = new double[N]; 
		
		for (int n=0; n< N; n++) {
	
		    Fr[n] = n*del_r*f[n];
	
		}
	
		//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		// Make auxiliary vector, Fr2 = [0,Fr,0,reverse(-Fr)], to evaluate DST with DFT
		//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		
		double[] Fr2 = new double[2*(N+1)];
	
		Fr2[0]=0; 
		Fr2[N+1] = 0;
		for (int i = 0; i<N; i++) {
			Fr2[i+1] = Fr[i];
		    Fr2[2*N+1-i] = -Fr[i];    
		}
	
		//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		// Perform DFT on auxilliary vector Fr2
		//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
		double[] Fk2 = new double[2*(N+1)];
		
		double[] Fk2i = new double[2*(N+1)]; // should be all zeros
		
		Fk2 = Fr2; 
		
		FastFourierTransform fourier = new FastFourierTransform();
		
		FastFourierTransform.BACKUP=false;
		fourier.setData(Fk2, Fk2i);
		fourier.transform();
		Fk2 = fourier.getImaginary();
	
		//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		//% Extract Fk from Fk2 = -1/N*[#,Fk,#,#]
		//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
		double[] Fk = new double[N];
		
		for (int i = 0; i<N; i++) {
			
			Fk[i] = -(N)*Fk2[i+1]; 
			
		}
		
		//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		//% Compute fk from Fk
		//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
	
		double [] fk = new double [N+1];
	
		// Fk(0) (k=0) does NOT equal zero, but we cannot calculate Fk(0)/0.
	
		
		for (int i = 1; i<N+1; i++) {
		  
		    
		    // Special consideration of zeroth mode, fk(k=0):
		     
		    fk[0] = fk[0] + 4.0*Math.PI*(Fr[i-1]*(i-1)*del_r)*del_r;
		    
		    fk[i] = 4.0*Math.PI*(Fk[i-1]/((i)*del_k))*del_r; 
	
		}
		
		return fk;
	}
	
	public double[] reverse(double[] fk, double del_r, double del_k) {
		

		//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		//% Make auxiliary vector, Fk = k*fk  
		//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

		int N = fk.length - 1;
		double[] Fk = new double[N];
		
		for (int i=0;i<N;i++) {
		    
		    Fk[i] = fk[i+1]*(i+1)*del_k; 
		    
		}

		//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		// Make auxiliary vector, Fk2 = [0,Fk,0,reverse(-Fk)], to evaluate DST with DFT
		//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		
		double[] Fk2 = new double[2*(N+1)];

		Fk2[0] = 0;
		Fk2[N+1] = 0; 
		
		for (int i=0; i<N; i++) {
			Fk2[i+1] = Fk[i]; 
			Fk2[2*N+1-i] = -Fk[i];    
		}
		    
		//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		//% Perform DFT on auxiliary vector Fk2
		//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
		double[] Fr2 = new double[2*(N+1)];
		
		double[] Fk2b = new double[2*(N+1)]; // should be all zeros
		
		FastFourierTransform fourier = new FastFourierTransform();
		
		FastFourierTransform.BACKUP=false;
		fourier.setData(Fk2, Fk2b);
		fourier.transform();
		Fr2 = fourier.getImaginary();

		//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		//% Extract Fr from Fr2
		//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		
		double[] Fr = new double[N];
		
		for (int i=0; i<N; i++) {

			Fr[i] = -(N)*Fr2[i+1]; 
		}

		//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		//% Compute fr from Fr=r*fr;
		//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

		double[] fr = new double[N];

		for (int i=1; i<N; i++) {
		   
		    //Special consideration for zeroth mode
		    fr[0] = fr[0] + 1/(2*Math.PI*Math.PI)*(Fk[i]*i*del_k)*del_k;

		    fr[i] = 1/(2*Math.PI*Math.PI)*(Fr[i]/(i*del_r))*del_k;
		    
		}
		
		return fr;
	}

	

}
