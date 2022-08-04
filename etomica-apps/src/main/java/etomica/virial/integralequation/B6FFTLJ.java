/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial.integralequation;

import etomica.math.numerical.SineTransform;
import etomica.potential.P2LennardJones;
import etomica.space.Space;
import etomica.space3d.Space3D;

/**
 * 
 * Modification of FFT Program  of B4 
 * @Srihari
 *
 */

public class B6FFTLJ {
	
public static void main(String[] args) {
	
double[] temps = {0.625,1.0,1.3,1.4,1.5,2.5,5.0,10.0};
	
for (int t=0; t<temps.length; t++) {
	
	double reducedTemp = temps[t]; // kT/epsilon
	
	int power = 14; // Defines discretization
	double r_max = 250; // Defines range of separation distance, r = [0 rmax]
	int N = (int) Math.pow(2, power);  // Number of grid points in the discretizations of r- and k-space
	
        double del_r = r_max/(N-1);
        double del_k = Math.PI/r_max;
		
		double[] fr = getfr( N, del_r,reducedTemp);
		
		SineTransform dst = new SineTransform();

		double[] fk = dst.forward(fr, del_r);
		
		double[] ffk 	= new double[N];
		double[] fffk	= new double[N];
		
		double[] ffffk  = new double[N];
		double[] fffffk = new double[N];
		
		for (int i = 0;i<N; i++) {
			ffk[i]    = fk[i]*fk[i];
			fffk[i]   = fk[i]*fk[i]*fk[i];
			ffffk[i]  = fk[i]*fk[i]*fk[i]*fk[i];
			fffffk[i] = fk[i]*fk[i]*fk[i]*fk[i]*fk[i];
		}
		
		double[] d1r = dst.reverse(fffk, del_r);

		double[] ffr    = dst.reverse(ffk, del_r);
		double[] fffr   = dst.reverse(fffk, del_r);
		
		double[] d2b1r = new double[N];
		for (int i = 1;i<(N-1); i++) {
			
			d2b1r[i] = fr[i]*ffr[i];
		}
		
		double[] d2b1k = dst.forward(d2b1r, del_r);
		
		double[] d2bk = new double[N];
		double[] f9k = new double[N];
		
		for (int i = 0;i<N; i++) {
			
			d2bk[i] = fk[i] * d2b1k[i];
			f9k[i] = d2b1k[i] * d2b1k[i];
		}
		
		double[] d2br = dst.reverse(d2bk, del_r);
		double[] f9r = dst.reverse(f9k, del_r);		
		
		double C3 = 0;
		double D1 = 0;
		double D2A = 0;
		double D2B = 0;
		
		double F1 = 0;
		double F2 = 0;
		double F3 = 0;
		double F4 = 0;
		double F5 = 0;
		double F6 = 0;
		double F7 = 0;
		double F8 = 0;
		double F9 = 0;
		double F10 = 0;
		double F11 = 0;
		double F12 = 0;
		double F13 = 0;
		double F14 = 0;
		double F15 = 0;
				
		for (int i = 1;i<(N-1); i++) {
			
			C3 = C3 + fr[i]*ffr[i]*(i*del_r)*(i*del_r)*del_r;
			 
			D1 = D1 + fr[i]*d1r[i]*(i*del_r)*(i*del_r)*del_r;
			
			D2A = D2A + fr[i]*(ffr[i]*ffr[i])*(i*del_r)*(i*del_r)*del_r;
			
			D2B = D2B + fr[i]*(d2br[i] )*(i*del_r)*(i*del_r)*del_r;
			
	  	}
		
		for (int i = 1;i<(N-1); i++){
			
			F1 = F1 + fffr[i] * fffr[i]*(i*del_r)*(i*del_r)*del_r;
			                         
			F2 = F2 + fffr[i] * ffr[i] * ffr[i]*(i*del_r)*(i*del_r)*del_r;
			
			F4 = F4 + ffr[i] * ffr[i] * fffr[i]* fr[i]*(i*del_r)*(i*del_r)*del_r;
			
			F3 = F3 + ffr[i] * ffr[i] * ffr[i] * ffr[i]*(i*del_r)*(i*del_r)*del_r;
			
			F5 = F5 + ffr[i] * ffr[i] * ffr[i] * ffr[i]*(i*del_r)*(i*del_r)*del_r*fr[i];
			
			F6 = F6 + fr[i] * fffr[i] * fffr[i]*(i*del_r)*(i*del_r)*del_r;
			
			F7 = F7 + d2br[i] * d2br[i]*(i*del_r)*(i*del_r)*del_r;
			
			F8 = F8 + d2br[i] * fffr[i]*(i*del_r)*(i*del_r)*del_r;
			
			F9 = F9 + ffr[i] * f9r[i]*(i*del_r)*(i*del_r)*del_r;  //f9r  = ((f*f)f) * ((f*f)f)  
			
			F10 = F10 +  d2br[i] * d2br[i] * fr[i]*(i*del_r)*(i*del_r)*del_r;
			
			F11 = F11 + ffr[i] * ffr[i] * d2br[i]*(i*del_r)*(i*del_r)*del_r;
			
			F12 = F12 +  d2br[i] * fr[i] * ffr[i] * ffr[i] *(i*del_r)*(i*del_r)*del_r;
			
			F13 = F13 + fr[i] * ffr[i] * f9r[i]*(i*del_r)*(i*del_r)*del_r;
			
			F14 = F14 +  d2br[i] * d2br[i] * fr[i]*(i*del_r)*(i*del_r)*del_r;
			
			F15 = F15 + fr[i] * d2br[i] * fffr[i]*(i*del_r)*(i*del_r)*del_r;  

		}//if(temps[t]==1) System.out.println(F3+ "  " +F5);
		
		C3 = C3 + 0.5*fr[N-1]*ffr[N-1]*(N-1)*del_r*(N-1)*del_r*del_r;
		
		D1 = D1 + 0.5*fr[N-1]*d1r[N-1]*(N-1)*del_r*(N-1)*del_r*del_r; 
		
		D2A = D2A + 0.5*fr[N-1]*( ffr[N-1]*ffr[N-1] )*(N-1)*del_r*(N-1)*del_r*del_r; 
		
		D2B = D2B + 0.5*fr[N-1]*( d2br[N-1] )*(N-1)*del_r*(N-1)*del_r*del_r; 
		
		F1 = F1 + 0.5*fffr[N-1] * fffr[N-1]*(N-1)*(del_r)*(N-1)*(del_r)*del_r;
        
		F2 = F2 + 0.5*fffr[N-1] * ffr[N-1] * ffr[N-1]*(N-1)*(del_r)*(N-1)*(del_r)*del_r;
		
		F4 = F4 + 0.5*ffr[N-1] * ffr[N-1] * fffr[N-1]* fr[N-1]*(N-1)*(del_r)*(N-1)*(del_r)*del_r;
		
		F3 = F3 + 0.5*ffr[N-1] * ffr[N-1] * ffr[N-1] * ffr[N-1]*(N-1)*(del_r)*(N-1)*(del_r)*del_r;
		
		F5 = F5 + 0.5*ffr[N-1] * ffr[N-1] * ffr[N-1] * ffr[N-1]*(N-1)*(del_r)*(N-1)*(del_r)*del_r*fr[N-1];
		
		F6 = F6 + 0.5*fr[N-1] * fffr[N-1] * fffr[N-1]*(N-1)*(del_r)*(N-1)*(del_r)*del_r;
		
		F7 = F7 + 0.5*d2br[N-1] * d2br[N-1]*(N-1)*(del_r)*(N-1)*(del_r)*del_r;
		
		F8 = F8 + 0.5*d2br[N-1] * fffr[N-1]*(N-1)*(del_r)*(N-1)*(del_r)*del_r;
		
		F9 = F9 + 0.5*ffr[N-1] * f9r[N-1]*(N-1)*(del_r)*(N-1)*(del_r)*del_r;  //f9r  = ((f*f)f) * ((f*f)f)  
		
		F10 = F10 +  0.5*d2br[N-1] * d2br[N-1] * fr[N-1]*(N-1)*(del_r)*(N-1)*(del_r)*del_r;
		
		F11 = F11 + 0.5*ffr[N-1] * ffr[N-1] * d2br[N-1]*(N-1)*(del_r)*(N-1)*(del_r)*del_r;
		
		F12 = F12 +  0.5*d2br[N-1] * fr[N-1] * ffr[N-1] * ffr[N-1] *(N-1)*(del_r)*(N-1)*(del_r)*del_r;
		
		F13 = F13 +0.5* fr[N-1] * ffr[N-1] * f9r[N-1]*(N-1)*(del_r)*(N-1)*(del_r)*del_r;
		
		F14 = F14 + 0.5* d2br[N-1] * d2br[N-1] * fr[N-1]*(N-1)*(del_r)*(N-1)*(del_r)*del_r;
		
		F15 = F15 + 0.5*fr[N-1] * d2br[N-1] * fffr[N-1]*(N-1)*(del_r)*(N-1)*(del_r)*del_r;  
		
		//if(temps[t]==1) System.out.println(0.5*ffr[N-1] * ffr[N-1] * ffr[N-1] * ffr[N-1]*(N-1)*(del_r)*(N-1)*(del_r)*del_r + " (---)  " +0.5*ffr[N-1] * ffr[N-1] * ffr[N-1] * ffr[N-1]*(N-1)*(del_r)*(N-1)*(del_r)*del_r* fr[N-1]);
		
		C3 = -1.0/(3.0)*4.0*Math.PI*C3;
		D1 = -3.0/(8.0)*4.0*Math.PI*D1;
		D2A = -3.0/(8.0)*4.0*Math.PI*D2A;
		D2B = -3.0/(8.0)*4.0*Math.PI*D2B;
		
		F1  = -5.0/(12.0)*4.0*Math.PI*F1;
		F2  = -5.0/(4.0)*4.0*Math.PI*F2;
		F3  = -5.0/(48.0)*4.0*Math.PI*F3;
		F4  = -5.0/(4.0)*4.0*Math.PI*F4;
		F5  = -5.0/(48.0)*4.0*Math.PI*F5;
		F6  = -5.0/(4.0)*4.0*Math.PI*F6;
		F7  = -5.0/(4.0)*4.0*Math.PI*F7;
		F8  = -5.0/(2.0)*4.0*Math.PI*F8;
		F9  = -5.0/(2.0)*4.0*Math.PI*F9;
		F10 = -5.0/(2.0)*4.0*Math.PI*F10;
		F11 = -5.0/(2.0)*4.0*Math.PI*F11;
		F12 = -5.0/(2.0)*4.0*Math.PI*F12;
		F13 = -5.0/(6.0)*4.0*Math.PI*F13;
		F14 = -5.0/(2.0)*4.0*Math.PI*F14;
		F15 = -5.0*4.0*Math.PI*F15;
				
		//double b = Standard.B2HS(1.0);
		
		/*double F[] = {F1,F2,F3,F4,F5,F6,F7,F8,F9,F10,F11,F12,F13,F14,F15};
		
		for(int i = 0 ;i<15;i++)
		{
			//System.out.println("F"+ (i+1) + "   " + F[i] );
			//System.out.println("   " + a);
				
		}*/
			
		System.out.println(reducedTemp + "   " + (F1+F2+F3+F4+F5+F6+F7+F8+F9+F10+F11+F12+F13+F14+F15));
			
		
	}
	}

	public static double[] getfr(int N, double del_r, double reducedTemp) {
		
		// Lennard-Jones Potential
	    double sigma = 1.0;
		double epsilon = 1.0;
		Space space = Space3D.getInstance();
		P2LennardJones p2 = new P2LennardJones(sigma, epsilon);
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
