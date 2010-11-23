package etomica.virial;

import etomica.potential.P2LennardJones;
import etomica.space.Space;
import etomica.space3d.Space3D;
import etomica.util.SineTransform;
import etomica.virial.cluster.Standard;

/**
 * 
 * Computes "simple" part of the f-bond-only formulation of B5.
 * 
 * @author Kate Shaul
 * @Modified for B5 - author Srihari
 *
 */


public class B5FFTLJ {
	
public static void main(String[] args) {
	
	double [] temps = new double[] { 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 2, 2.5, 3, 5, 10, 15, 20, 50, 100};
	
	for (int t=0; t<temps.length; t++) {
		
		double reducedTemp = temps[t]; // kT/epsilon
		
		int power = 14; // Defines discretization
		double r_max = 200; // Defines range of separation distance, r = [0 rmax]
		int N = (int) Math.pow(2, power);  // Number of grid points in the discretizations of r- and k-space
	
        double del_r = r_max/(N-1);
        double del_k = Math.PI/r_max;
		
		double[] fr = getfr( N, del_r,reducedTemp);
		
		SineTransform dst = new SineTransform();

		double[] fk = dst.forward(fr, del_r, r_max);
		
		double[] ffk = new double[N];
		double[] fffk = new double[N];
		double[] ffffk = new double[N];
		double[] fck = new double[N];
		
		for (int i = 0;i<N; i++) {
			ffk[i] = fk[i]*fk[i];
			fffk[i] = fk[i]*fk[i]*fk[i];
			ffffk[i] = fk[i]*fk[i]*fk[i]*fk[i];
		}
		
		double[] ffffr = dst.reverse(ffffk, del_r, r_max);
		double[] fffr  = dst.reverse(fffk, del_r, r_max);
		double[] ffr   = dst.reverse(ffk, del_r, r_max);
		
		double[] E6ar = new double[N];
		double[] E6br = new double[N];
		double[] E7gr = new double[N];
		double[] E7acr = new double[N];				
		
		for (int i = 1;i<(N-1); i++) {
			
			E6ar[i] = ffr[i]*fr[i]*fffr[i];
			E6br[i] = ffr[i]*ffr[i]*ffr[i];
			E7gr[i] = ffr[i]*ffr[i]*ffr[i]*fr[i];
			E7acr[i] = fr[i]*ffr[i];
			
		}
		
		double[] E7ack = dst.forward(E7acr, del_r, r_max);
		double[] E7ac1k = new double[N];
		
		for (int i = 0;i<N; i++) {
			
			E7ac1k[i] = fk[i]*E7ack[i];
		}
		
		double[] E7ac1r = dst.reverse(E7ac1k, del_r, r_max);
		
		for (int i = 1;i<(N-1); i++) {
			
				E7ac1r[i] = fr[i]*ffr[i]*E7ac1r[i];
		}
		
		double E5 = 0;
		double E6B = 0;
		double E6A = 0;
		double E7A =0;
		double E7G = 0;
		
		for (int i = 1;i<(N-1); i++) {
						 
			E5 = E5 + fffr[i]*ffr[i]*(i*del_r)*(i*del_r)*del_r;
			
			E6B = E6B + E6br[i]*(i*del_r)*(i*del_r)*del_r;
		
			E6A = E6A + E6ar[i]*(i*del_r)*(i*del_r)*del_r;
			
			E7A = E7A + E7ac1r[i]*(i*del_r)*(i*del_r)*del_r;
			
			E7G = E7G + E7gr[i]*(i*del_r)*(i*del_r)*del_r;
		}
		
	    	E5 = E5 + 0.5*( fffr[N-1]*ffr[N-1] )*(N-1)*del_r*(N-1)*del_r*del_r; 
	    	
	    	E6B = E6B + 0.5*( E6br[N-1] )*(N-1)*del_r*(N-1)*del_r*del_r;
	    	
	    	E6A = E6A + 0.5*( E6ar[N-1] )*(N-1)*del_r*(N-1)*del_r*del_r;
	    	
	    	E7A = E7A + 0.5*( E7ac1r[N-1] )*(N-1)*del_r*(N-1)*del_r*del_r;
	    	
	    	E7G = E7G + 0.5*( E7gr[N-1] )*(N-1)*del_r*(N-1)*del_r*del_r;

		
		E5 = -2.0/(5.0)*4.0*Math.PI*E5;	
		
		E6A = -2.0*4.0*Math.PI*E6A;
		
		E6B = -1.0/(3.0)*4.0*Math.PI*E6B;
		
		E7A = -2.0*4.0*Math.PI*E7A;
		
		E7G = -1.0/(3.0)*4.0*Math.PI*E7G;
		
		double b = Standard.B2HS(1.0);

	/*	System.out.println("B3  (T* = "+ reducedTemp + ") = "+ B3);
		System.out.println("D1* (T* = "+ reducedTemp + ") = "+ D1/(b*b*b));
		System.out.println("D2A*(T* = "+ reducedTemp + ") = "+ D2A/(b*b*b));
		System.out.println("D2B*(T* = "+ reducedTemp + ") = "+ D2B/(b*b*b));
		System.out.println("D2* (T* = "+ reducedTemp + ") = "+ (D2A+D2B)/(b*b*b));
		
		System.out.println(reducedTemp + "   " + D1/(b*b*b) + "    " + (D2A+D2B)/(b*b*b));*/
		
	/*	System.out.println(reducedTemp + "   " + E5/(b*b*b*b));
		System.out.println(reducedTemp + "   " + E6A/(b*b*b*b));
		System.out.println(reducedTemp + "   " + E6B/(b*b*b*b));
		System.out.println(reducedTemp + "   " + E7A/(b*b*b*b));
		System.out.println(reducedTemp + "   " + E7G/(b*b*b*b));*/
		
		//System.out.println(reducedTemp + "   " + E5/(b*b*b*b) + "    " + E6A/(b*b*b*b)+ "    " + E6B/(b*b*b*b)+ "    " + E7A/(b*b*b*b)+ "    " + E7G/(b*b*b*b));
		
		//System.out.println(reducedTemp + "   " + (E5+E6A+E6B+E7A+E7G));
		
		double[] D = {E5,E6A,E6B,E7A,E7G};
		
		//if(temps[t]==1)
		{
			for(int i = 0 ;i<5;i++)
			{
				//System.out.println("D" + (i+1) + "   " +D[i]);
			}
			
			System.out.println(temps[t] + "   " + (E5+E6A+E6B+E7A+E7G));
		}
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
