/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial.integralequation;

import etomica.potential.P2LennardJones;
import etomica.space.Space;
import etomica.space3d.Space3D;

/**
 * 
 * Combines the Percus-Yevick and hypernetted-chain approximations to the virial coefficients to 
 * generate a self-consistent (thermodynamically consistent) approximation.
 * 
 * @author kate
 *
 */


public class BnSCALJ {
	
public static void main(String[] args) {
        
		// To make sure that everything is working fine:
		System.out.println("Literature values for LJ with sigma = 1 and T* = 1");
		System.out.println("B2 = -5.3158 (Sun & Teja 1996)");
		System.out.println("B3 =  1.8849 (Sun & Teja 1996) ");
		System.out.println("B4PY(c) = -2.9394 (Dyer et al 2001)");
		System.out.println("B4PY(c) = -2.925 (imprecise value from Henderson, Kim, and Oden (1966))");
		System.out.println("B4PY(v) = -6.627 (imprecise value from Henderson, Kim, and Oden (1966))");
		System.out.println("B4HNC(c) = -11.4784 (Dyer et al 2001) ");
		System.out.println("B4HNC(v) = -4.39 (Dyer et al 2001) ");
		System.out.println("B4  = -2.463(7) (MSMC)");
		System.out.println("B5PY(c) = -37.402 (Barker et al 1966)\n");
		
		
		
		int power = 15; // Defines discretization
		
		double r_max = 100; // Defines range of separation distance, r = [0 rmax]
		
		int N = 1<<power;
		double del_r = r_max/(N-1);
		
		int m = 4; // highest order of virial coefficient to be calculated
		double reducedTemp;
		if (args.length == 0) {
		}
		else if (args.length == 4) {
			power = Integer.parseInt(args[0]);
			r_max = Double.parseDouble(args[1]);
			reducedTemp = Double.parseDouble(args[2]);
			m = Integer.parseInt(args[3]);
		} else {
        	throw new IllegalArgumentException("Incorrect number of arguments passed.");
        }
		
		
		double [] reducedTemps = new double[] { 0.6, 0.8, 1.0, 1.05, 1.1, 1.15, 1.2, 1.3, 1.4, 1.6, 2, 2.5, 3, 5, 10, 15, 20, 30, 50, 100, 500};
		//double [] reducedTemps = new double[] { 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 2, 2.5, 3, 5, 10, 15, 20, 50, 100};

		//double [] reducedTemps = new double[] { 1.0};
		
		System.out.println("temp   phi    B"+m+"SCA:" );
		for (int i=0;i<reducedTemps.length;i++) {
            reducedTemp = reducedTemps[i];

            // Number of grid points in the discretizations of r- and k-space


            // Get Mayer function for this discretization
            double[] fr = getfr(N, del_r, reducedTemp);
            double[] rdfdr = getrdfdr(N, del_r, reducedTemp);

            PercusYevick pyc = new PercusYevick();
            pyc.setRoute(true);
            pyc.setrdfdr(rdfdr);

            double[] BPYC = pyc.computeB(fr, m, N, del_r, false);


            HypernettedChain hncc = new HypernettedChain();

            hncc.setRoute(true);
            hncc.setrdfdr(rdfdr);

            double[] BHNCC = hncc.computeB(fr, m, N, del_r, false);

            PercusYevick pyv = new PercusYevick();
            pyv.setRoute(false);
            pyv.setrdfdr(rdfdr);
			
			double[] BPYV = pyv.computeB(fr, m, N, del_r, false);
			
			HypernettedChain hncv = new HypernettedChain(); 
			
			hncv.setRoute(false);
			hncv.setrdfdr(rdfdr);
			
			double[] BHNCV = hncv.computeB(fr, m, N, del_r, false);
			
			double BmPYC = BPYC[m-2];
			double BmHNCC = BHNCC[m-2];
			double BmPYV = BPYV[m-2];
			double BmHNCV = BHNCV[m-2];
			
			double delPY = BmPYC - BmPYV;
			double delHNC = BmHNCC - BmHNCV;
			
			double phi = delPY/(delPY-delHNC);
			
			double BmSCA = (1.0-phi)*BmPYC + phi*BmHNCC;
			
			//double BmSCA = (1.0-phi)*BmPYV + phi*BmHNCV;
			
			
			
			System.out.println(reducedTemp + "    " + phi + "   " + BmSCA );
			
			
			
			//System.out.println(reducedTemp + "  " + B[5-2]);
		}
		
		

	}

	public static double[] getfr(int N, double del_r, double reducedTemp) {
	    
		double sigma = 1.0;
		
	    double epsilon = 1.0;
		
		Space space = Space3D.getInstance();
		
		P2LennardJones p2 = new P2LennardJones(sigma, epsilon);
		
		
		
		double[] fr = new double[N];  // Holds discretization of Mayer function in r-space
		
		for (int n = 0; n<N; n++) {
			
			double r = n*del_r;
			
			double u = p2.u(r*r);
			
			double x = -u/reducedTemp;
			
			if ( Math.abs(x) < 0.01) {
				
				fr[n] = x + x*x/2.0 + x*x*x/6.0 + x*x*x*x/24.0 + x*x*x*x*x/120.0;	
				
			} else {
				
				fr[n] = Math.exp(x)-1.0;
				
			}
	
		}
		
		return fr;
		
	}
	
	public static double[] getrdfdr(int N, double del_r, double reducedTemp) {
	    
		double sigma = 1.0;
		
	    double epsilon = 1.0;
		
		Space space = Space3D.getInstance();
		
		P2LennardJones p2 = new P2LennardJones(sigma, epsilon);
		
		
		
		double[] rdfdr = new double[N];  // Holds discretization of Mayer function in r-space
		
		for (int n = 1; n<N; n++) {
			
			double r = n*del_r;
			
			double u = p2.u(r*r);
			
			double x = -u/reducedTemp;
			
			//To match Henderson and Oden (1966):
			//rdfdr[n] = (-Math.exp(x)/reducedTemp)*p2.du(r*r);
			// My gn and hn include the factor of e; above would double count it:
			rdfdr[n] = (-1.0/reducedTemp)*p2.du(r*r);

		}
		
		return rdfdr;
		
	}


}
