package etomica.virial;

import etomica.potential.P2ArgonAziz1993;
import etomica.potential.P2LennardJones;
import etomica.potential.Potential2SoftSpherical;
import etomica.space.Space;
import etomica.space3d.Space3D;
import etomica.util.SineTransform;


/**
 * 
 * This method uses iteration and B3SphericallySymmetric.java to compute B3 for a spherically symmetric potential to within a given tolerance by 
 * either FFT or FFT and one layer of quadrature.
 * 
 * @author kate
 *
 */


public class B3ForSphericallySymmetricU {
	
	public B3ForSphericallySymmetricU() {
		
	}
	
public static void main(String[] args) {

	boolean printapalooza = false;
	
	boolean allFFT = false;

	int power = 6;
	
	double r_max = 100;
	
	Space space = Space3D.getInstance();
	
	boolean argon = false;
	
	double[] temps;
	
	Potential2SoftSpherical p2; 
	
	if (argon) {
		p2 = new P2ArgonAziz1993(space);
		//P2QChemInterpolated p2 = new P2QChemInterpolated(space);
	
		// Temperatures for argon potentials:
		temps = new double[] { 100, 133.15, 150, 200, 250, 300, 400, 500, 600, 700, 800, 900, 1000 }; // Kelvin
		// Some of the temperatures considered by Mas, Lotrich, and Szalewicz (1999):  
		//double[] temps = { 113.15, 133.15, 150.65, 203.15, 323.15, 423.15, 573.16, 673.16, 773.15, 923.15 }; // Kelvin
		// Primarily temperatures considered by Malejevsky et al:  
		//double[] temps = { 130, 135, 140,145,150,155, 160, 170, 180, 190, 200, 220, 250, 265, 280, 295, 310, 325, 340, 398, 423, 473, 573, 673};

	} else {
		
		p2 = new P2LennardJones(space,1.0,1.0);
		temps =  new double[] { 1.0 }; 
	}
	
	if (allFFT) {
		System.out.println("For B3, FFT employed to compute integral over r12\n");
	} else {
		System.out.println("For B3, Trapezoidal rule employed to compute integral over r12\n");
	}
	
	System.out.println();
	System.out.println("T (K)    B2      abs(B2(r_max)-B2(r_max/2))     B3       abs(B3(N)-B3(N/2))     log2(N)     rmax");
	System.out.println();
	
	for (int t=0; t<temps.length; t++) {
		
		double temp = temps[t];
   
		double[] results = getConvergence((Potential2SoftSpherical) p2, power, r_max, temp,  printapalooza, allFFT);

		System.out.println(temp + "    "  + results[0] + "    " + results[1] + "    "+ results[3] + "    " + results[4] + "   " + results[5] + "   " + results[2]);
	
	}
	

	}

	public static double[] computeB(double[] fr, int N, double r_max, boolean allFFT) {
		
		/*******************************************
		/*******************************************
		 * The method relies upon Fourier transforms.  The 3-D Fourier transform simplifies to the following form upon asserting spherical symmetry:
		 * 
		 * f(k) = integral over (constant*f(r)*r*sin(k*r)) dr
		 * 
		 * Following Duh and Haymet (1995), the transform becomes an an even simpler sine transform if one uses auxiliary functions of the form, F(r) = r*f(r)
		 * 
		 * f(k) = integral over (constant*F(r)*sin(k*r))dr
		 * 
		 * We do not use F(k) = k*f(k) = (k*constant*F(r)*sin(k*r))dr because the k=0 values are very important to us: We need them to be nonzero and defined to calculate the virial coefficients.
		 * 
		 * We use a fast transform method to achieve high accuracy for LJ potentials at low T.
		 * 
		/*******************************************
		/********************************************/
	
		double[] B = new double[2];
		
		double del_r = r_max/(N-1);

		SineTransform dst = new SineTransform();
		
		double[] dummy = fr;
		
		double[] fk = dst.forward(dummy, del_r, r_max);
		
		double B2 = -1.0/(2.0)*(fk[0]);
		
		B[0] = B2;
			
		double[] ffk = new double[N];
		for (int i = 0;i<N; i++) {
			ffk[i] = fk[i]*fk[i];
		}
			
		dummy = ffk;
		double[] ffr = dst.reverse(dummy, del_r, r_max);
			
		// Trapezoidal rule:
		
		double B3 = 0;
		
		for (int i = 1;i<(N-1); i++) {
		 
			B3 = B3 + fr[i]*ffr[i]*(i*del_r)*(i*del_r)*del_r;
		
		}
			
		B3 = B3 + 0.5*fr[N-1]*ffr[N-1]*r_max*r_max*del_r; 
		
		B3 = -1.0/(3.0)*4.0*Math.PI*B3;
			
		B[1] = B3;
	
		return B;
	}

	public static double[] getConvergence (Potential2SoftSpherical p2, int power, double r_max, double temp, boolean printapalooza, boolean allFFT) {
		
	
		double [] results = new double[6];
		double[] newB;
		
		double[] B = new double [2];
		double[] errors = new double [2];
		double error = 1.0;
		boolean molPerL = false;
		
		/* 
		 * Two variables need to be explored: rmax and the number of grid points, N.
		 * 
		 */
		
		int powerNew = power;
		double r_maxNew = r_max;
		int N = (int) Math.pow(2, power);
		
		double tol = 1e-3;
		
		while (error > tol) {  // increase r_max to improve accuracy of B2
			
			r_max = r_maxNew;
			
			powerNew = power;
			B[0] = 0;
			B[1] = 0;
			
			while (error > tol) {  // increase N to improve accuracy of B3
				
				powerNew = powerNew+1;
				
				N = (int) Math.pow(2, powerNew);
				
				double del_r = r_max/((double)(N-1));
	
				double[] fr = getfr( p2, N, del_r,temp);	
				
				newB = computeB(fr, N, r_max, allFFT);
				
				if (molPerL) {
					newB[0]=newB[0]*0.60221415;
					newB[1]=newB[1]*0.60221415*0.60221415;
				}
				
				error = Math.abs(B[1]-newB[1]);
				
				errors[1] = error;
				
				B = newB;
				
				//System.out.println(power2 + "   " + B[0] +  "  " + error);
				
				if (printapalooza) {
					
					System.out.println(temp + "    "  + B[0] + "    " + errors[0] + "    "+ B[1] + "    " + errors[1] + "   " + powerNew + "   " + r_max);
				}
				
				if (powerNew > 16) {
					break;
				}
				
			}
			
			if (r_max < 300) {
				r_maxNew = r_max*2;
			}
			
			int NNew = N*2;
		
			double del_r = r_maxNew/((double)(NNew-1));
	
			double[] fr = getfr(p2, NNew, del_r,temp);
	
			newB = computeB(fr, NNew, r_maxNew, allFFT);
			
			if (molPerL) {
				newB[0]=newB[0]*0.60221415;
				newB[1]=newB[1]*0.60221415*0.60221415;
			}			
			
			error = Math.abs(B[0]-newB[0]);
			
			
			errors[0] = error;
		
			B = newB;
			
			if (printapalooza) {
				
				System.out.println(temp + "    "  + B[0] + "    " + errors[0] + "    "+ B[1] + "    " + errors[1] + "   " + powerNew + "   " + r_max);
			}
				
		}
		
		results[0] = B[0];
		results[1] = errors[0];
		results[2] = r_maxNew;
		results[3] = B[1];
		results[4] = errors[1];
		results[5] = powerNew;
		
		return results;
	}

	public static double[] getfr(Potential2SoftSpherical p2, int N, double del_r, double temp) {
		
		double[] fr = new double[N];  // Holds discretization of Mayer function in r-space
		
		for (int n = 0; n<N; n++) {
			
			double r = n*del_r;
			
			// The following if statement only makes sense for argon potentials:
			// if ( r < 2.5) {
				
				//fr[n]=-1.0;
				
			//} else {
				
				double u = p2.u(r*r); // returned in Kelvin

				double x = -u/temp;
				
				if ( Math.abs(x) < 0.01) {
					
					fr[n] = x + x*x/2.0 + x*x*x/6.0 + x*x*x*x/24.0 + x*x*x*x*x/120.0;
					
					
				} else {
					
					fr[n] = Math.exp(x)-1.0;
					
				}
				
			///}
			
		}
		
		return fr;
		
	}

}
