/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial.integralequation;

import etomica.potential.*;
import etomica.space.Space;
import etomica.space3d.Space3D;
import etomica.units.Kelvin;
import etomica.util.ParameterBase;
import etomica.util.ParseArgs;

/**
 * 
 * This is just a main method for PercusYevick.java.  It creates a
 * discretization of the He Mayer function, which PercusYevick then employs to
 * compute Percus-Yevick approximations of the virial coefficients up to mth
 * order.  The calculation will increase the number of discretization points
 * and maximum separation until the coefficient converges to within the desired
 * tolerance.
 * 
 * @author kate
 * @author Andrew Schultz
 */
public class BnPYCHe {

    enum PotentialChoice {
        SIMPLE, PCKLJS, PCJS
    }

    public static class BnPYHeParams extends ParameterBase {
        public int n = 2;
        public double T = 500;
        public boolean verbose = false;
        public double tol = 1e-8;
        public PotentialChoice potentialChoice = PotentialChoice.PCJS;
        public double rmax = 40;
        public double sigma = 0;
        public boolean classical = true;
    }
    
    public static void main(String[] args) {
        BnPYHeParams params = new BnPYHeParams();
        boolean isCommandline = args.length > 0;
        if (isCommandline) {
            ParseArgs.doParseArgs(params, args);
        }
        else {
            params.n = 2;
            params.T = 100;
            params.verbose = false;
            params.potentialChoice = PotentialChoice.PCJS;
            params.rmax = 40;
            params.tol = 1e-8;
            params.sigma = 0;
        }
        int m = params.n;
        double T = Kelvin.UNIT.toSim(params.T);
        double tol = params.tol;
        boolean verbose = params.verbose;
        double sigma = params.sigma;
        Space space = Space3D.getInstance();
        if (!params.classical && params.potentialChoice == PotentialChoice.PCJS) {
            throw new RuntimeException("Don't know how to do semiclassical PCJS");
        }
        PotentialChoice pc = params.potentialChoice;
        P2HeSimplified p2Simple = new P2HeSimplified(space);
        P2HePCKLJS p2PCKLJS = new P2HePCKLJS(space, sigma);
        P2HePCJS p2PCJS = new P2HePCJS(space, sigma);
        Potential2SoftSpherical p2 = pc == PotentialChoice.SIMPLE ? p2Simple : (pc == PotentialChoice.PCKLJS ? p2PCKLJS : p2PCJS);
        Potential2Spherical p2sc = pc == PotentialChoice.SIMPLE ? new P2HeSimplified(space).makeQFH(T) : new P2HePCKLJS(space, sigma).makeQFH(T);
        double rMax = params.rmax;
        while (true) {
            double[][] results = getConvergence(params.classical ? p2 : p2sc, m, rMax, tol, T, verbose, pc == PotentialChoice.SIMPLE ? 1.6 : 0);
            if (results[2][0] > tol) {
                double newRMax = results[2][2];
                if (newRMax > 2.5*rMax) {
                    if (tol > 2*params.tol) {
                        tol = 2*params.tol;
                    }
                    rMax *= 2.0;
                    continue;
                }
                tol *= 2;
                if (verbose) System.out.println("tol => "+tol);
                continue;
            }
            if (params.verbose) {
                for (int i=0; i<m-1; i++) {
                    System.out.println("B"+(i+2)+"  "+results[0][i]+" "+results[1][i]);
                }
            }
            else {
                System.out.println(String.format("%20.15e  %10.4e  %2d  %4d", results[0][m-2], results[1][m-2], (int)results[2][1], (int)results[2][2]));
            }
            break;
        }
    }

    public static double[][] getConvergence (Potential2Spherical p2, int m, double r_max, double tol, double kT, boolean printapalooza, double core) {
        
        double[][] results = new double[3][m];
        if (m==2) {
            results = new double[3][3];
        }
        double[] errorBPY = new double[m];
        double[] oldBPY = new double[m];
        double[] powerErr = new double[m];
        double[] B = null;
    
        int N = 0;
        
        int power = 8;
        double error = 2.0*tol;
        boolean searchForPower = true;
        while (error > tol) {  // increase N to improve accuracy of B2 and B3
            
            N = (int) Math.pow(2, power);
            if (m == 2) {
                N++;
                double del_r = r_max/(N-1);
                if (core > 0 && del_r > core/2)  {
                    power++;
                    continue;
                }
                double myRMax = r_max;
                if (core > 0) {
                    int iCore = (int)(core / del_r);
                    del_r *= core / (iCore*del_r);
                    myRMax = del_r*(N-1);
                }
                double[] fr = getfr(p2, N, del_r, kT, core);
                double sum = 0;
                for (int i = 1; i<N-1; i++) {
                    double r = i*del_r;
                    sum += r*r*fr[i];
                }
                sum += 0.5*myRMax*myRMax*fr[N-1];
                B = new double[]{-2*Math.PI*sum*del_r};
            }
            else {
            
                double del_r = r_max/(N-1);
                if (core > 0 && del_r > core/2)  {
                    power++;
                    continue;
                }
                if (core > 0) {
                    int iCore = (int)(core / del_r);
                    del_r *= core / (iCore*del_r);
                }

                PercusYevick py = new PercusYevick();
                py.setRoute(true);
                
                double[] fr = getfr(p2, N, del_r, kT, core);
                B = py.computeB(fr, m, N, del_r, false);
            }

            error=0;
            for (int i = 0; i<m-1; i++) {
                errorBPY[i] = Math.abs(B[i]-oldBPY[i]);
            }
            error = errorBPY[m-2];
    
            
            if (printapalooza) {
                System.out.print(power+" "+r_max);
                for (int i = 0; i<m-1; i++) {   
                    System.out.print("    "  + B[i] + "    " + errorBPY[i]);
                }
                System.out.println();
            }

            if (error <= tol && searchForPower) {
                searchForPower = false;
                error = 1;
                System.arraycopy(errorBPY, 0, powerErr, 0, m-1);
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

            for (int i = 0; i<m-1; i++) {   
                oldBPY[i] = B[i];
            }
                
        }
        
        for (int i = 0; i<m-1; i++) {   
            results[0][i] = B[i];
            results[1][i] = Math.sqrt(errorBPY[i]*errorBPY[i]+powerErr[i]*powerErr[i]);
        }
        
        results[2][0] = error;
        results[2][1] = power;
        results[2][2] = r_max;
        
        return results;
    }

    public static double[] getfr(Potential2Spherical p2, int N, double del_r, double temp, double core) {
        
        double[] fr = new double[N];  // Holds discretization of Mayer function in r-space
        
        for (int n = 0; n<N; n++) {
            
            double r = n*del_r;
            if (r > core-del_r*0.1 && r < core+del_r*0.1) {
                r = core;
            }
            double u = p2.u(r*r);
            double x = -u/temp;
            
            if ( Math.abs(x) < 0.01) {
                fr[n] = x*(1.0 + x*(1.0 + x*(1.0 + x*(1.0 + x*(1.0 + x/6.0)/5.0)/4.0)/3.0)/2.0);
            }
            else {
                fr[n] = Math.exp(x)-1.0;
            }
            if (r == core) {
                fr[n] = 0.5*(-1 + fr[n]);
            }
        }

        return fr;
    }
}
