package etomica.virial;

/**
 * 
 * This is just a main method for PercusYevick.java...
 * 
 * It creates a discretization of the hard-sphere Mayer function, which PercusYevick then employs to compute Percus-Yevick
 * approximations of the virial coefficients up to (m+1) order.  The second and third coefficients are fully accurate (the 
 * PY approximation is exact).
 * 
 * @author kate
 *
 */


public class BnPYHardSphere {
	
public static void main(String[] args) {
        
		System.out.println("Perera computed the following values for sigma = 1");
		System.out.println("B2 = 2.0943951");
		System.out.println("B3 = 2.74155678");
		System.out.println("B4 = 2.72740397");
		System.out.println("B5 = 2.33000141");
		System.out.println("B6 = 1.81030163");
		System.out.println("B7 = 1.31877804");
		System.out.println("B8 = 0.91708435");
		System.out.println("B9 = 0.61576568");
		System.out.println("B10 = 0.40227821");
		System.out.println("B11 = 0.2570954");
		System.out.println("B12 = 0.16137562\n");
		
		int power = 18; // Defines discretization
		double r_max = 25; // Defines range of separation distance, r = [0 rmax]
		double sigma = 1; // Defines hard-sphere system
		
		if (args.length == 0) {
		}
		else if (args.length == 2) {
			power = Integer.parseInt(args[0]);
			r_max = Double.parseDouble(args[1]);
		} else {
        	throw new IllegalArgumentException("Incorrect number of arguments passed.");
        }
		
		// Number of grid points in the discretizations of r- and k-space
		int N = (int) Math.pow(2, power) - 1;
    	
		// m+1 is the highest order of virial coefficient to be calculated
        int m = 5; 
        
        double del_r = r_max/(N-1);

		double r = 0.0;
	
		double[] fr = new double[N];  // Holds discretization of Mayer function in r-space
		
		for (int n = 0; n<N; n++) {

			if (r <= sigma) {
				fr[n] = -1.0;
			}  else {
				fr[n] = 0.0;
			}
				
			r += del_r; 

		}
		
		PercusYevick py = new PercusYevick(); 
		double[] B = py.computeB(fr, m, N, del_r);
		
		System.out.println("Values computed here:");
		for (int i=0;i<m;i++) {
			System.out.println("B" + (i+2) + " = " + B[i]);
		}
		

	}

}
