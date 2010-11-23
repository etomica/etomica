package etomica.virial;

import etomica.potential.P2LennardJones;
import etomica.space.Space;
import etomica.space3d.Space3D;
import etomica.util.SineTransform;
import etomica.virial.cluster.Standard;

/**
 * 
 * Computes "simple" part of the f-bond-only formulation of B4.
 * 
 * @author Kate Shaul
 *
 */


public class B4FFTLJ {
	
public static void main(String[] args) {
	
	double [] temps = new double[] { 0.6, 0.8, 1.0, 1.05, 1.1, 1.15, 1.2, 1.3, 1.4, 1.6, 2, 2.5, 3, 5, 10, 15, 20, 30, 50, 100, 500};
	
	for (int t=0; t<temps.length; t++) {
		
		double reducedTemp = temps[t]; // kT/epsilon
		
		int power = 14; // Defines discretization
		double r_max = 100; // Defines range of separation distance, r = [0 rmax]
		int N = (int) Math.pow(2, power);  // Number of grid points in the discretizations of r- and k-space
	
        double del_r = r_max/(N-1);
        double del_k = Math.PI/r_max;
		
		double[] fr = getfr( N, del_r,reducedTemp);
		
		SineTransform dst = new SineTransform();

		double[] fk = dst.forward(fr, del_r, r_max);
		
		double[] ffk = new double[N];
		double[] fffk = new double[N];
		for (int i = 0;i<N; i++) {
			ffk[i] = fk[i]*fk[i];
			fffk[i] = fk[i]*fk[i]*fk[i];
		}
		
		double[] d1r = dst.reverse(fffk, del_r, r_max);

		double[] ffr = dst.reverse(ffk, del_r, r_max);
		
		double[] d2b1r = new double[N];
		for (int i = 1;i<(N-1); i++) {
			
			d2b1r[i] = fr[i]*ffr[i];
		}
		
		double[] d2b1k = dst.forward(d2b1r, del_r, r_max);
		
		double[] d2bk = new double[N];
		for (int i = 0;i<N; i++) {
			
			d2bk[i] = fk[i]*d2b1k[i];
		}
		
		double[] d2br = dst.reverse(d2bk, del_r, r_max);
		
		double B3 = 0;
		double D1 = 0;
		double D2A = 0;
		double D2B = 0;
		for (int i = 1;i<(N-1); i++) {
			
			B3 = B3 + fr[i]*ffr[i]*(i*del_r)*(i*del_r)*del_r;
			 
			D1 = D1 + fr[i]*d1r[i]*(i*del_r)*(i*del_r)*del_r;
			
			D2A = D2A + fr[i]*(ffr[i]*ffr[i])*(i*del_r)*(i*del_r)*del_r;
			D2B = D2B + fr[i]*(d2br[i] )*(i*del_r)*(i*del_r)*del_r;
			
		}
		
	    B3 = B3 + 0.5*fr[N-1]*ffr[N-1]*(N-1)*del_r*(N-1)*del_r*del_r;
		
		D1 = D1 + 0.5*fr[N-1]*d1r[N-1]*(N-1)*del_r*(N-1)*del_r*del_r; 
		
		D2A = D2A + 0.5*fr[N-1]*( ffr[N-1]*ffr[N-1] )*(N-1)*del_r*(N-1)*del_r*del_r; 
		
		D2B = D2B + 0.5*fr[N-1]*( d2br[N-1] )*(N-1)*del_r*(N-1)*del_r*del_r; 
		

		B3 = -1.0/(3.0)*4.0*Math.PI*B3;
		
		D1 = -3.0/(8.0)*4.0*Math.PI*D1;
		
		D2A = -3.0/(8.0)*4.0*Math.PI*D2A;
		
		D2B = -3.0/(8.0)*4.0*Math.PI*D2B;
		
		double b = Standard.B2HS(1.0);

		//System.out.println("B3  (T* = "+ reducedTemp + ") = "+ B3);
		//System.out.println("D1* (T* = "+ reducedTemp + ") = "+ D1/(b*b*b));
		//System.out.println("D2A*(T* = "+ reducedTemp + ") = "+ D2A/(b*b*b));
		//System.out.println("D2B*(T* = "+ reducedTemp + ") = "+ D2B/(b*b*b));
		//System.out.println("D2* (T* = "+ reducedTemp + ") = "+ (D2A+D2B)/(b*b*b));
		
		System.out.println(reducedTemp + "   " + (D1 + (D2A+D2B)) );
		
	}
	}

	public static double[] getfr(int N, double del_r, double reducedTemp) {
		
		// Lennard-Jones Potential
	    double sigma = 1.0;
		double epsilon = 1.0;
		Space space = Space3D.getInstance();
		P2LennardJones p2 = new P2LennardJones(space, sigma, epsilon);
		double u;
		
		double r = 0.0;
		
		double[] fr = new double[N];  // Holds discretization of Mayer function in r-space
		
		for (int n = 0; n<N; n++) {
			
			u = p2.u(r*r);
			
			double x = -u/reducedTemp;
			
			if ( Math.abs(x) < 0.01) {
				
				fr[n] = x + x*x/2.0 + x*x*x/6.0 + x*x*x*x/24.0 + x*x*x*x*x/120.0;
				
				
			} else {
				
				fr[n] = Math.exp(x)-1.0;
				
			}
			
			// fr[n] = Math.exp(x)-1.0;
			
			
			r += del_r; 
	
		}
		
		return fr;
		
	}


}
