/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial;

import etomica.math.numerical.SineTransform;
import etomica.potential.IPotential2;

/**
 * This is a main method to compute B3, dB3/dT and d2B3/dT2 for a spherical potential via FFT
 * 
 * It creates a discretization of the Mayer function, and computes B3 via FFT.
 * The tolerance can be set.  Initial values of n and rMax will be increased to find convergence
 * to the desired tolerance.  The function may not converge if n and rMax are too far away from
 * appropriate values.
 */
public class B3FFT {

    /**
     * @param p2            the pair potential
     * @param n             the number of points used for the discretization is 2^n
     * @param rMax          the truncation distance of the potential
     * @param tol0          the desired tolerance for the virial coefficient
     * @param kT            the temperature
     * @param verbose       if true, any increase in the tolerance will be noted
     * @param printapalooza if true, results from each rMax and n will be printed
     * @return              returns the results
     *                      results[0] are the B3 and its derivatives
     *                      results[1] are the estimates of error in B3 and its derivatives
     *                      results[2] are the error in B3, the final n and final rMax used
     */
    public static double[][] go(IPotential2 p2, int n, double rMax, double tol0, double kT, boolean verbose, boolean printapalooza) {
        double tol = tol0;
        while (true) {
            double[][] results = getConvergence(p2, n, rMax, tol, kT, printapalooza);
            if (results[2][0] > tol) {
                double newRMax = results[2][2];
                if (newRMax > 2.5*rMax) {
                    if (tol > 2*tol0) {
                        tol = 2*tol0;
                    }
                    rMax *= 2.0;
                    continue;
                }
                tol *= 2;
                if (verbose) System.out.println("tol => "+tol);
                continue;
            }
            return results;
        }
    }

    public static double[][] getConvergence(IPotential2 p2, int n, double r_max, double tol, double kT, boolean printapalooza) {
        int d = 0;

        double[] errorB = new double[3];
        double[] oldB = new double[3];
        double[] powerErr = new double[3];
        double[] B = null;

        int power = 8;
        double error = 2.0 * tol;
        boolean searchForPower = true;
        while (error > tol) {  // increase N to improve accuracy of B2 and B3
            int N = (int) Math.pow(2, power);
            double del_r = r_max/(N-1);

            // Get Mayer function (f) and df/dT for this discretization
            double[] fr = getfr(p2, N, del_r, kT);
            double[] dfdT = getdfdT(p2, N, del_r, kT);
            double[] d2fdT2 = getd2fdT2(p2, N, del_r, kT);
//        for (int i=0; i<fr.length; i++) {
//            System.out.println(i*del_r+" "+fr[i]+" "+dfdT[i]+" "+d2fdT2[i]);
//        }
//        System.exit(0);

            if (n == 2) {
                double sum = 0, dsum = 0, d2sum = 0;
                for (int i = 1; i<N; i++) {
                    double r = i*del_r;
                    sum += r*r*fr[i];
                    dsum += r*r*dfdT[i];
                    d2sum += r*r*d2fdT2[i];
                }
                sum *= -2*Math.PI*del_r;
                dsum *= -2*Math.PI*del_r;
                d2sum *= -2*Math.PI*del_r;
                B = new double[]{sum, dsum, d2sum};
            }
            else {
                double[] fk = SineTransform.forward(fr, del_r);

                double[] fk2 = new double[fk.length];
                for (int i = 0; i < fk.length; i++) {
                    fk2[i] = fk[i] * fk[i];
                }
                double[] ffr = SineTransform.reverse(fk2, del_r);
                double d2sum = 0;
                fk = SineTransform.forward(dfdT, del_r);

                for (int i = 0; i < fk.length; i++) {
                    fk2[i] = fk[i] * fk[i];
                }
                double[] ffrdT = SineTransform.reverse(fk2, del_r);
                for (int i = 0; i < fr.length; i++) {
                    d2sum += (ffrdT[i] * fr[i] * i) * i;
                }
                d2sum *= 2;
                // dB/dT is the cluster integral with one one bond as df/dT and a coefficient
                // of -1 instead of -1/3
                double sum = 0, dsum = 0;
                for (int i = 0; i < fr.length; i++) {
                    sum += (ffr[i] * fr[i] * i) * i;
                    dsum += (ffr[i] * dfdT[i] * i) * i;
                    d2sum += (ffr[i] * d2fdT2[i] * i) * i;
                }
                double fourpidr3 = 4.0 * Math.PI * (del_r * del_r * del_r);
                sum *= -fourpidr3 / 3.0;
                dsum *= -fourpidr3;
                d2sum *= -fourpidr3;
                B = new double[]{sum, dsum, d2sum};
            }


            for (int i = 0; i<B.length; i++) {
                errorB[i] = Math.abs(B[i]-oldB[i]);
            }
            error = errorB[d];

            if (printapalooza) {
                System.out.print(power+" "+r_max);
                for (int i = 0; i<B.length; i++) {
                    System.out.print("    "  + B[i] + "    " + errorB[i]);
                }
                System.out.println();
            }

            if (error <= tol && searchForPower) {
                searchForPower = false;
                error = 1;
                System.arraycopy(errorB, 0, powerErr, 0, errorB.length);
            }

            if (error > tol) {
                if (power > 20) {
                    if (printapalooza) {
                        System.out.println("power > 17, perhaps you should try lowering the tolerance");
                    }
                    break;
                }
                power++;
                if (!searchForPower) {
                    //still double number of points, but also extend range
                    r_max *= 2;
                }
            }

            System.arraycopy(B, 0, oldB, 0, B.length);

        }

        double[][] results = new double[3][3];

        for (int i = 0; i<B.length; i++) {
            results[0][i] = B[i];
            results[1][i] = Math.sqrt(errorB[i]*errorB[i]+powerErr[i]*powerErr[i]);
        }

        results[2][0] = error;
        results[2][1] = power;
        results[2][2] = r_max;

        return results;
    }

    /**
     * Returns a discretized f function containing N points with a
     * separation interval of del_r
     */
    public static double[] getfr(IPotential2 p2, int N, double del_r, double reducedTemp) {

		double[] fr = new double[N];  // Holds discretization of Mayer function in r-space
		
		for (int n = 0; n<N; n++) {
			double r = n*del_r;
			double u = p2.u(r*r);
			double x = -u/reducedTemp;
			if ( Math.abs(x) < 0.01) {
                fr[n] = x + x * x / 2.0 + x * x * x / 6.0 + x * x * x * x / 24.0 + x * x * x * x * x / 120.0;
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
    public static double[] getdfdT(IPotential2 p2, int N, double del_r, double reducedTemp) {
        double[] dfdT = new double[N];  // Holds discretization of Mayer function in r-space

        // eqn below is inf*0 for r=0
        dfdT[0] = 0;

        for (int n = 1; n < N; n++) {
            double r = n * del_r;
            double u = p2.u(r * r);
            double x = -u / reducedTemp;
            dfdT[n] = Math.exp(x);
            if (dfdT[n] != 0) {
                dfdT[n] *= -x / reducedTemp;
            }
        }

        return dfdT;
    }

    /**
     * Returns a discretized df/dT function containing N points with a
     * separation interval of del_r
     */
    public static double[] getd2fdT2(IPotential2 p2, int N, double del_r, double reducedTemp) {
        double[] d2fdT2 = new double[N];  // Holds discretization of Mayer function in r-space

        // eqn below is inf*0 for r=0
        d2fdT2[0] = 0;

        for (int n = 1; n < N; n++) {
            double r = n * del_r;
            double u = p2.u(r * r);
            double x = -u / reducedTemp;
            d2fdT2[n] = Math.exp(x);
            if (d2fdT2[n] != 0) {
                d2fdT2[n] *= (x * x + 2 * x) / (reducedTemp * reducedTemp);
            }
            if (Double.isNaN(d2fdT2[n])) {
                throw new RuntimeException("oops "+r+" "+u+" "+x);
            }
        }

        return d2fdT2;
    }
}
