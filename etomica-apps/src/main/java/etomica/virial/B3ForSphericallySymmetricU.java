/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial;

import etomica.math.numerical.SineTransform;
import etomica.potential.IPotential2;
import etomica.potential.P2QChemInterpolated;
import etomica.space.Space;
import etomica.space3d.Space3D;
import etomica.util.ParameterBase;


/**
 * 
 * Compute B3 for a spherically symmetric potential to within a given tolerance, specified with getConvergece(), by 
 * either FFT or FFT and one layer of quadrature.
 * 
 * @author kate
 *
 */



public class B3ForSphericallySymmetricU {
	
	public B3ForSphericallySymmetricU() {
		
	}
	
public static void main(String[] args) {

	
	DampingParams params = new DampingParams();
    
   
	if (args.length == 6 ) {
		params.a1 = Integer.parseInt(args[0]);
		params.a2 = Integer.parseInt(args[1]);
		params.Rvdw = Double.parseDouble(args[2]);
		params.basis = Integer.parseInt(args[3]);
		params.fixedRvdw = Boolean.parseBoolean(args[4]);
		params.tempSet = Integer.parseInt(args[5]);
    } 
	
	int a1 = params.a1;
	int a2 = params.a2;
	double Rvdw = params.Rvdw;
	int basis = params.basis;
	boolean fixedRvdw = params.fixedRvdw;
	int tempSet = params.tempSet;
	
	boolean printapalooza = false;
	
	boolean allFFT = false;
	
	boolean qm = true;

	// Minimum log2(N) and maximum separation distance to be considered:
	int power = 10;
	double r_max = 200;
	
	Space space = Space3D.getInstance();
	
	boolean argon = true;
	
	double[] temps;
	
	
	P2QChemInterpolated p2 = new P2QChemInterpolated();
	p2.setDampingParams(a1,a2,Rvdw,basis, fixedRvdw);
	p2.setDisp(true);
	p2.setSCF(true);
	p2.initialize();
	
	
	//P2ArgonTangAndToennies2003 p2 = new P2ArgonTangAndToennies2003(space);
	//Potential2SoftSpherical p2 = new P2ArgonAziz1993(space);
	
	if (argon) {
		//
	   
		//
		//System.out.println("BJ w /aug-cc-pV"+basis+"Z, a1 = " + a1 + ", a2 = " + a2 + "\n");
		// Temperatures for argon potentials:
		if (tempSet == 2) {
			temps = new double[] { 100, 133.15, 150, 200, 250, 300, 400, 500, 600, 700, 800, 900, 1000 }; // Kelvin
		} else {
			// Primarily temperatures considered by Malejevsky et al:  
			temps = new double[] { 130, 135, 140,145,150,155, 160, 170, 180, 190, 200, 220, 250, 265, 280, 295, 310, 325, 340, 398, 423, 473, 573, 673};
		}
		
		// Some of the temperatures considered by Mas, Lotrich, and Szalewicz (1999):  
		//double[] temps = { 113.15, 133.15, 150.65, 203.15, 323.15, 423.15, 573.16, 673.16, 773.15, 923.15 }; // Kelvin
		

	} else {
		
		//p2 = new P2LennardJones(space,1.0,1.0);
		temps =  new double[] { 1.0 }; 
	}
	
	if (allFFT) {
		System.out.println("For B3, FFT employed to compute integral over r12\n");
	} else {
		System.out.println("For B3, Trapezoidal rule employed to compute integral over r12\n");
	}
	
	
	
	System.out.println();
	if (qm) {
		System.out.println("T(K)    B2      B2C      abs(B2C(N)-B2C(N/2))     	B2QM	B3       abs(B3(N)-B3(N/2))     log2(N)     rmax	del_r");
	} else {
		System.out.println("T(K)    B2      abs(B2(N)-B2(N/2))     B3       abs(B3(N)-B3(N/2))     log2(N)     rmax");
	}
	
	System.out.println();
	
	for (int t=0; t<temps.length; t++) {
		
		double temp = temps[t];
   
		double[] results = getConvergence(p2, power, r_max, temp,  printapalooza, allFFT);

		double B2 = results[0];;
		double B2Error = results[1];
		double r_maxF = results[2];
		double B3 = results[3];
		double B3Error = results[4];
		double powerF = results[5];
		double NF = Math.pow(2,powerF);
		double del_rF = r_maxF/(NF-1);
		
		if (qm) {
			double qmB2 = computeB2QM(p2, temp);
			double totalB2 = qmB2 + B2;
			System.out.println(temp + "    "  + totalB2 + "    " + B2 + "    "+ B2Error + "    "+ qmB2 +"    " +B3 + "    " + B3Error + "   " + powerF + "   " + r_maxF+ "   " +del_rF);
			//System.out.println(temp + "    "+ qmB2 );
		} else {
			
			
			System.out.println(temp + "    "  + B2 + "    " + B2Error + "    "+ B3 + "    " + B3Error + "   " + powerF + "   " + r_maxF);
		}
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
		
		double[] fk = dst.forward(dummy, del_r);
		
		double B2 = -1.0/(2.0)*(fk[0]);
		
		B[0] = B2;
			
		double[] ffk = new double[N];
		for (int i = 0;i<N; i++) {
			ffk[i] = fk[i]*fk[i];
		}
			
		dummy = ffk;
		double[] ffr = dst.reverse(dummy, del_r);
			
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
	
	public static double computeB2QM(IPotential2 p2, double temp) {

		
		 double del_r = 0.0001;
		 double h = 6.62606896e-34; //J*s
		 double k = 1.3806504e-23; // J/K
	     double m = 39.948*(1.660538782e-27); //kg
	     double constant =  h*h/(24*Math.PI*m*k*k*k*temp*temp*temp);
	     // constant [=] J*J*s*s/(kg*J*J*J) [=] s*s/(J*kg) 
	     // 	J = kg*m*m/(s*s); (s*s) = kg*m*m/J
	     // constant [=] (kg*m*m/J)/(J*kg) = m*m/(J*J)
	     // but u will be given in Kelvin and r in Angstroms
	     
	     constant = constant*(1e20)*k*k; // A*A/(K*K)

	     double r12=2.5+del_r;
	     double u12Backward=p2.u(r12*r12);
	     double u12;
	     double slope;
	     double e12;
	     double qmB2 = 0;
	     while(r12 < 10) {
	    	 
	    	 r12=r12+del_r;
			 
			 u12 = p2.u(r12*r12);
			 
			 slope = (u12-u12Backward)/del_r; // K/A

			 e12 = Math.exp(-u12/temp);

		     qmB2 = qmB2 + e12*slope*slope*r12*r12*del_r; // K*K/(A*A)*A*A*A
		     
	    	 u12Backward=u12;
	    	 
	    	
	    	
	     }

	     
	     u12 = p2.u(r12*r12);
	     slope = (u12-u12Backward)/del_r;
	     e12 = Math.exp(-u12/temp);
	     qmB2 = qmB2 + 0.5*e12*slope*slope*r12*r12*del_r;
	     
	     //System.out.println(qmB2*constant*0.60221214);
	     //System.exit(0);

	     return qmB2*constant*0.60221214;  // from A^3/molecule to cm^3/mole, 6.022e23/(1e24)

	}

	public static double[] getConvergence (IPotential2 p2, int power, double r_max, double temp, boolean printapalooza, boolean allFFT) {
		
		r_max = 200;

		double [] results = new double[6];
		double [] newB;
		double B2Old=0;double B2 = 1;double B2Error=1;
		double B3Old=0;double B3 = 1;double B3Error=1;

		int N = (int) Math.pow(2, power);

		double error = 1.0;
		double tol = 0.01;
		while (error > tol) {  // increase N to improve accuracy of B2 and B3
			
			N = (int) Math.pow(2, power);
			
			double del_r = r_max/((double)(N-1));

			double[] fr = getfr( p2, N, del_r,temp);	
			
			newB = computeB(fr, N, r_max, allFFT);	
			B2 = newB[0]*0.60221415;
			B3 = newB[1]*0.60221415*0.60221415;
			
			B2Error = Math.abs(B2-B2Old);
			B3Error = Math.abs(B3-B3Old);
				
			error = Math.max(B2Error, B3Error);
				
				
			if (printapalooza) {
				
				System.out.println(temp + "    "  + B2 + "    " + B2Error + "    "+ B3 + "    " + B3Error + "   " + power + "   " + r_max);
			}
				
			power = power+1;
			B2Old=B2;
			B3Old=B3;
				
		}
		
		results[0] = B2;
		results[1] = B2Error;
		results[2] = r_max;
		results[3] = B3;
		results[4] = B3Error;
		results[5] = power;
		
		return results;
	}

	public static double[] getfr(IPotential2 p2, int N, double del_r, double temp) {
		
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
				
				/*
				if (Double.isNaN(fr[n])){
					System.out.println(r+" "+u*k/JPerHartree*1e6+" "+fr[n]);
				}
				*/
				
				
			///}
			
		}

		
		return fr;
		
	}
	
	 public static class DampingParams extends ParameterBase {
		    //TZ
		 /*
		 	public int a1 = 79;	        
	        public int a2 = 136;   
	        private double Rvdw = 3.65;
	        private int basis = 3;
	        private boolean fixedRvdw = false;
		 	*/
	        
		 	//DZ
	    	public int a1 = 80;	        
	        public int a2 = 149;   
	        private double Rvdw = 3.688;
	        private int basis = 2;
	        private boolean fixedRvdw = false;
	        
	    
	        
	        public int tempSet = 2; 
	  }
	 
	 protected static double k = 1.3806503e-23; //J/K   1.3806503e-23
	    protected static double JPerHartree = 4.359744e-18;  //4.359744e-18

}
