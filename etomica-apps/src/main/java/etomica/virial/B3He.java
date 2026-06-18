/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial;

import etomica.potential.P2HePCJS;
import etomica.units.Kelvin;
import etomica.util.Constants;
import etomica.util.ParameterBase;
import etomica.util.ParseArgs;

/**
 * This is a main method to compute B3 and dB3/dT for Lennard-Jones via FFT
 * 
 * It creates a discretization of the Lennard-Jones Mayer function, and computes
 * B3 via FFT.  The resulting values are accurate to ~13 digits (from T=0.39 up to >> 1000)
 * 
 * @author Andrew Schultz
 */
public class B3He {

    public static class B3HeParams extends ParameterBase {
        public double T = 500;
        public boolean verbose = false;
        public double tol = 1e-8;
        public double rmax = 40;
        public double sigma = 0;
        public int n = 3;
    }

    public static void main(String[] args) {
        B3HeParams params = new B3HeParams();
        boolean isCommandline = args.length > 0;
        if (isCommandline) {
            ParseArgs.doParseArgs(params, args);
        }
        else {
            params.T = 100;
            params.verbose = false;
            params.rmax = 40;
            params.tol = 1e-8;
            params.sigma = 0;
            params.n = 3;
        }

        P2HePCJS p2 = new P2HePCJS();
        double TK1 = 10, TK2 = 16;
        int nT = 61;
        double kB = Constants.BOLTZMANN_K;
        double tol = params.tol;
        boolean verbose = params.verbose;
        int n = params.n;

        for (int iT=0; iT<nT; iT++) {
            double TK = TK1 + (TK2-TK1) / (nT - 1) * iT;
            double reducedTemp = Kelvin.UNIT.toSim(TK);
            double rMax = params.rmax;

            double[][] results = B3FFT.go(p2, n, rMax, tol, reducedTemp, verbose, false);
            if (params.verbose) {
                for (int i=0; i<results[0].length; i++) {
                    System.out.println("B3d"+(i)+"  "+results[0][i]+" "+results[1][i]);
                }
            }
            else {
                System.out.printf("%4.1f ", TK);
                for (int i = 0; i < results[0].length; i++) {
                    System.out.print(" " + results[0][i]*Math.pow(kB,i));
                }
                System.out.printf("  %2d  %4d\n", (int)results[2][1], (int)results[2][2]);
            }

        }
    }
}
