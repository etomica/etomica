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
		
		int power = 17; // Defines discretization
		double r_max = 100; // Defines range of separation distance, r = [0 rmax]
		int N = 1<<power;  // Number of grid points in the discretizations of r- and k-space
	
        double del_r = r_max/(N-1);
		
		double[] fr = getfr( N, del_r,reducedTemp);
		
		SineTransform dst = new SineTransform();

		double[] fk = dst.forward(fr, del_r, r_max);
		
		double[] ffk = new double[N];
		double[] fffk = new double[N];
		for (int i = 0;i<N; i++) {
			ffk[i] = fk[i]*fk[i];
			fffk[i] = fk[i]*fk[i]*fk[i]; 
		}
		
		double[] ffr = dst.reverse(ffk, del_r, r_max);
		
		double[] ffk2 = dst.forward(ffr, del_r, r_max);
		
		double[] fffk2 = new double [N];
		for (int i = 0;i<N; i++) {
			fffk2[i] = ffk2[i]*fk[i]; 
		}
		
		/* Using fffk2 to compute d1r is more accurate than using fffk: the result converges *to the same answer* more quickly.
		 * This indicates that error in the transform is compounded in the fffk computation, and that some of this error is somehow 
		 * removed when taking the reverse and and forward transforms of ffk.  It seems likely 
		 * not all of the error is removed by this ah hoc procedure, which would indicate that even faster convergence is possible.
		 */
		
		double[] d1r = dst.reverse(fffk2, del_r, r_max);
		
		double[] d2b1r = new double[N];
		for (int i = 0;i<N; i++) {
			
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
		
		double[] cB3 = new double[N];
		double[] cB4Approx = new double[N];
		cB3[0] = fr[0]*ffr[0];
		cB4Approx[0] = fr[0]*( d1r[0]+ (ffr[0]*ffr[0]) + d2br[0] );
		
		
		for (int i = 1;i<(N-1); i++) {
			
			double c = (i*del_r)*(i*del_r)*del_r;
			
			B3 = B3 + fr[i]*ffr[i]*c;
			 
			D1 = D1 + fr[i]*d1r[i]*c;
			
			D2A = D2A + fr[i]*(ffr[i]*ffr[i])*c;
			
			D2B = D2B + fr[i]*(d2br[i] )*c;
			
			cB3[i] = fr[i]*ffr[i];
			    
			cB4Approx[i] = fr[i]*( d1r[i] + (ffr[i]*ffr[i]) + d2br[i] );
			
		}
		
		double c = r_max*r_max*del_r;
		
	    B3 = B3 + 0.5*fr[N-1]*ffr[N-1]*c;
		
		D1 = D1 + 0.5*fr[N-1]*d1r[N-1]*c; 
		
		D2A = D2A + 0.5*fr[N-1]*( ffr[N-1]*ffr[N-1] )*c; 
		
		D2B = D2B + 0.5*fr[N-1]*( d2br[N-1] )*c; 
		
		cB3[N-1] = fr[N-1]*ffr[N-1];
	    
		cB4Approx[N-1] = fr[N-1]*( d1r[N-1] + (ffr[N-1]*ffr[N-1]) + d2br[N-1]);
		
		
			
		

		B3 = -1.0/(3.0)*4.0*Math.PI*B3;
		
		D1 = -1.0/(8.0)*4.0*Math.PI*D1;
		
		D2A = -1.0/(8.0)*4.0*Math.PI*D2A;
		
		D2B = -1.0/(8.0)*4.0*Math.PI*D2B;
		
		
		double [] c1k = dst.forward(cB3, del_r, r_max);
		
		double [] c2k = dst.forward(cB4Approx, del_r, r_max);
			
		double B3F = -1.0/(3.0)*(c1k[0]); // B3 for m = 1
		
		double B4PY = -1.0/(4.0)*(c2k[0]); // B3 for m = 1
		
		//double b = Standard.B2HS(1.0);
		//System.out.println("B3  (T* = "+ reducedTemp + ") = "+ B3);
		//System.out.println("D1* (T* = "+ reducedTemp + ") = "+ D1/(b*b*b));
		//System.out.println("D2A*(T* = "+ reducedTemp + ") = "+ D2A/(b*b*b));
		//System.out.println("D2B*(T* = "+ reducedTemp + ") = "+ D2B/(b*b*b));
		//System.out.println("D2* (T* = "+ reducedTemp + ") = "+ (D2A+D2B)/(b*b*b));
		
		//System.out.println(reducedTemp + "   " + B3 );
		//System.out.println(reducedTemp + "   " + B3F );
		
		
		//B4-B4Bridge
		//System.out.println(reducedTemp + "   " + 3*(D1 + (D2A+D2B)) );
		
		//B4PY
		System.out.println(reducedTemp + "   " + 2*(D1 + (D2A+D2B)) );
		//System.out.println(reducedTemp + "   " + B4PY );
		
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
