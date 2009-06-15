package etomica.virial;

import etomica.potential.P2LennardJones;
import etomica.space.Space;
import etomica.space3d.Space3D;

public class B3LJ {
	
	public B3LJ() {
		
	}
	
public static void main(String[] args) {
	
		int power = 17;
		double reducedTemp = 1.0; // kT/epsilon
		double r_max = 10;
		boolean allFFT = true;
		
		if (args.length == 0) {
			
		} else if (args.length == 3) { 
			
			reducedTemp = Double.parseDouble(args[0]);
			power = Integer.parseInt(args[1]);
			r_max = Double.parseDouble(args[2]);
			
			
		} else if (args.length == 6) {
			reducedTemp = Double.parseDouble(args[0]);
			power = Integer.parseInt(args[1]);
			r_max = Double.parseDouble(args[2]);
			
			allFFT = Boolean.parseBoolean(args[5]);
			
		} else {
        	throw new IllegalArgumentException("Incorrect number of arguments passed to B3LJ.");
        }
		
		
		// Number of grid points in the discretizations of r-space;
        // There are N+1 grid points in k-space.
		int N = (int) Math.pow(2, power) - 1;
    	
		// M+1 is the highest order of virial coefficient to be calculated
        int M = 2; 
        
        
		
		double del_r = r_max/(N-1);

		double[] fr = getfr( N, del_r,reducedTemp);
		
		
		
		B3SphericallySymmetric b3 = new B3SphericallySymmetric();
		System.out.println("B3LJ at T* = " + reducedTemp);
		
		
		if (allFFT) {
			System.out.println("FFT employed to compute integral over r12\n");
		} else {
			System.out.println("Riemann sum employed to compute integral over r12\n");
		}
		
		double[] B = b3.computeB(fr, M, N, del_r, allFFT);
		System.out.println(power + "   " + r_max + "   " + del_r + "   " + B[0] + "  "+ B[1]);
		
		double[] newB;
		
		double error2 = 1.0;
		double r_max2 = r_max;
		while (error2 > 1e-5) {
			
			r_max2 = r_max2 + 2;
			del_r = r_max2/(N-1);
			fr = getfr( N, del_r,reducedTemp);
			
			newB = b3.computeB(fr, M, N, del_r, allFFT);
			
			//double diffB2 = Math.abs(B[0]-newB[0]);
			double diffB3 = Math.abs(B[1]-newB[1]);
			//error2 = Math.max(diffB2, diffB3);
			error2 = diffB3;
			
			B = newB;
			System.out.println(power + "   " + r_max2 + "   " + del_r + "   " + B[0] + "  "+ B[1]);
			
		}
		System.out.println("  \n");
		
		double error = 1.0;
		double power2 = power;
		while (error > 1e-6) {
			
			power2 = power2+1;
			N = (int) Math.pow(2, power2) - 1;
			fr = getfr( N, del_r,reducedTemp);
			
			newB = b3.computeB(fr, M, N, del_r, allFFT);
			
			//double diffB2 = Math.abs(B[0]-newB[0]);
			double diffB3 = Math.abs(B[1]-newB[1]);
			// error = Math.max(diffB2, diffB3);
			error = diffB3;
			
			B = newB;
			// System.out.println(power2 + "   " + r_max2 + "   " + del_r + "   " + B[0] + "  "+ B[1]);
			System.out.println(power2 + "   " + r_max2 + "   " + del_r + "   " + B[1]);
			
		}
		
		double maxPower = power2;
		double max_r_max = r_max2;
		if (args.length == 6) {
			maxPower = Integer.parseInt(args[3]);
			max_r_max = Double.parseDouble(args[4]);
		}
		
		System.out.println("log2(N): " + power +" to " + maxPower);
		
		System.out.println("r_max: " + r_max +" to " + max_r_max);
		
		
		double[] Bfinal = B;
		
		System.out.println("  \n");
		
		boolean convergencePlot = true;
		
		if (convergencePlot) {
			
			while (r_max <= max_r_max) {
				
				double powerP = power;
				N = (int) Math.pow(2, powerP) - 1;
				del_r = r_max/(N-1);
				fr = getfr( N, del_r,reducedTemp);
				
				newB = b3.computeB(fr, M, N, del_r, allFFT);
				
				// double diffB2 = Math.abs(newB[0]-Bfinal[0]);
				double diffB3 = Math.abs(newB[1]-Bfinal[1]);
				
				// System.out.print(r_max + "   " + del_r + "  "+ diffB2 + "  " + diffB3);
				System.out.print(r_max + "   " + del_r + "  "+ diffB3);
				
				
				while (powerP < maxPower) {
					
					powerP = powerP+1;
					N = (int) Math.pow(2, powerP) - 1;
					fr = getfr( N, del_r,reducedTemp);
					
					newB = b3.computeB(fr, M, N, del_r, allFFT);
					
					// diffB2 = Math.abs(newB[0]-Bfinal[0]);
					diffB3 = Math.abs(newB[1]-Bfinal[1]);
		
					System.out.print("  " + diffB3);
					//System.out.print("  "+ diffB2 + "  " + diffB3);
					
				}
				
				r_max = r_max + 2;
				
				System.out.println();
				
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
			
			fr[n] = Math.exp(-u/reducedTemp)-1.0;
			
			r += del_r; 

		}
		
		return fr;
		
	}

}
