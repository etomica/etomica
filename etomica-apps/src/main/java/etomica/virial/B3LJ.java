/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial;

import etomica.potential.P2LennardJones;
import etomica.space.Space;
import etomica.space3d.Space3D;
import etomica.math.numerical.SineTransform;

/**
 * This is a main method to compute B3 and dB3/dT for Lennard-Jones via FFT
 * 
 * It creates a discretization of the Lennard-Jones Mayer function, and computes
 * B3 via FFT.  The resulting values are accurate to ~13 digits (from T=0.39 up to >> 1000)
 * 
 * @author Andrew Schultz
 */
public class B3LJ {
	
    public static void main(String[] args) {
        double reducedTemp = 1.0;
        if (true) {
            if (args.length == 0) {
                System.err.println("Usage: B3LJ temp");
                System.exit(1);
            }
            reducedTemp = Double.parseDouble(args[0]);
        }

        double[] values = value(reducedTemp);
		System.out.println(values[0]+" "+values[1]);
	}
    
    public static double[] value(double temperature) {
        double r_max = 50; // Defines range of separation distance, r = [0 rmax)
        if (temperature>1000) r_max *= Math.pow(temperature/1000, -0.25);
        int power = 12;
        int N = 1<<power;
        double del_r = r_max/N;
        
        // Get Mayer function (f) and df/dT for this discretization
        double[] fr = getfr(N, del_r, temperature);
        double[] dfdT = getdfdT(N, del_r, temperature);
        
        SineTransform dst = new SineTransform();

        double[] fk = dst.forward(fr, del_r);

        double[] fk2 = new double[fk.length];
        for (int i=0; i<fk.length; i++) {
            fk2[i] = fk[i]*fk[i];
        }
        double[] ffr = dst.reverse(fk2, del_r);
        // dB/dT is the cluster integral with one one bond as df/dT and a coefficient
        // of -1 instead of -1/3
        double sum = 0, dsum = 0;
        for (int i=0; i<fk.length; i++) {
            sum += (ffr[i]*fr[i]*i)*i;
            dsum += (ffr[i]*dfdT[i]*i)*i;
        }
        double fourpidr3 = 4.0*Math.PI * (del_r*del_r*del_r);
        sum *= -fourpidr3/3.0;
        dsum *= -fourpidr3;
        return new double[]{sum,dsum};
    }

    /**
     * Returns a discretized f function containing N points with a
     * separation interval of del_r
     */
    public static double[] getfr(int N, double del_r, double reducedTemp) {
	    
		double sigma = 1.0;
		
	    double epsilon = 1.0;
		
		Space space = Space3D.getInstance();
		
		P2LennardJones p2 = new P2LennardJones(space, sigma, epsilon);
		
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


	/**
	 * Returns a discretized df/dT function containing N points with a
	 * separation interval of del_r
	 */
    public static double[] getdfdT(int N, double del_r, double reducedTemp) {
        double sigma = 1.0;
        double epsilon = 1.0;
        Space space = Space3D.getInstance();
        P2LennardJones p2 = new P2LennardJones(space, sigma, epsilon);
        double[] dfdr = new double[N];  // Holds discretization of Mayer function in r-space

        // eqn below is inf*0 for r=0
        dfdr[0] = 0;

        for (int n = 1; n<N; n++) {
            double r = n*del_r;
            double u = p2.u(r*r);
            double x = -u/reducedTemp;
            dfdr[n] = -x/reducedTemp*Math.exp(x);
        }
        
        return dfdr;
    }
}
