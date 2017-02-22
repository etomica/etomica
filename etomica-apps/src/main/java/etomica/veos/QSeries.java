/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.veos;

import etomica.math.SpecialFunctions;
import etomica.math.discrete.FixedSumIterator;

/**
 * Fluid properties based on direct expansion of the canonical partition function Q in terms
 * of irreducible cluster integrals beta.  Values of beta integrals are given by BetaSource
 * instance associated with this class.
 *
 */
public class QSeries {
    
    private BetaSource betaSource;
    
    public QSeries(BetaSource bs) {
        betaSource = bs;
    }

    /**
     * Computes the coefficient b_j based on the beta coefficients given by the BetaSource.
     */
    public double bCalc(int j) {
        int maxk = Math.min(betaSource.maxIndex(), j-1);
        FixedSumIterator iter = new FixedSumIterator(j-1);
        iter.setSum(j-1);
        iter.reset();
        int[] n;
        double b = 0;
        while((n = iter.next()) != null) {
            double prod = 1.0;
            for(int k=1; k<=maxk; k++) {
                if(n[k-1]==0) continue;//use k-1 because indexing starts at 0
                prod *= Math.pow(j*betaSource.beta(k), n[k-1])/SpecialFunctions.gamma(n[k-1]+1);
            }
            b += prod;
        }
        return b/(j*j);
    }
    
    /**
     * Computes coefficients of the polynomial in volume for Q for the given N.
     * @return An array qC, such that qC[k] is the coefficient of V^(k+1) for the Q series in volume
     */
    public double[] coeffQ(int N) {
        double[] b = new double[N];
        b[0] = 1.0;
        for(int j=2; j<=N; j++) {
            b[j-1] = bCalc(j);
        }
        FixedSumIterator iter = new FixedSumIterator(N);
        iter.setSum(N);
        iter.reset();
        int[] m;

        double[] q = new double[N];
        while((m = iter.next()) != null) {
            double prod = 1.0;
            int msum = 0;
            for(int j=1; j<=N; j++) {
                if(m[j-1]==0) continue;
                msum += m[j-1];
                prod *= Math.pow(b[j-1], m[j-1])/SpecialFunctions.gamma(m[j-1]+1);
            }
            q[msum-1] += prod;
        }
        return q;
    }
    
    /**
     * Computes Q * N!/V^N based on the volume series for Q given by qCalc.
     * Evaluates partition function for a set of densities given by the input rho, and returns Q/V^N 
     * for each density as the corresponding element in the return array.
     */
    public double[] qCalc(int N, double[] rho) {
        double[] qVals = new double[rho.length];
        double[] q = coeffQ(N);
        for(int i=0; i<rho.length; i++) {
            double V = N/rho[i];
            qVals[i] = 0;
            double vprod = Math.exp(SpecialFunctions.lnFactorial(N)-N*Math.log(V));
            for(int j=1; j<=q.length; j++) {
                vprod *= V;
                qVals[i] += q[j-1]*vprod;
            }
        }
        return qVals;
    }

    
    /**
     * Computes the pressure via dlnQ/dV based on the volume series for Q given by coeffQ.
     * Evaluates pressure for a set of densities given by the input rho, and returns the pressure 
     * for each density as the corresponding element in the return array.
     * @param T temperature, given simply to yield the pressure as T * dlnQ/dV
     */
    public double[] pCalc(int N, double[] rho, double T) {
        double[] p = new double[rho.length];
        double[] q = coeffQ(N);
        for(int i=0; i<rho.length; i++) {
            double V = N/rho[i];
            double dsum = 0, nsum = 0;
            double vprod = 1.0;
            for(int j=1; j<=q.length; j++) {
                nsum += j*q[j-1]*vprod;
                vprod *= V;
                dsum += q[j-1]*vprod;
            }
            p[i] = T*nsum/dsum;
        }
        return p;
    }
        
    public static void main(String[] args) {
        double T = 0.8 * 8.0 / 27.0;
        QSeries ush = new QSeries(new VanderWaals(T,20));
        int N = 10;
        for(int j=2; j<=10; j++) {
            System.out.println(j+" "+ush.bCalc(j));
        }
        System.out.println();
        double[] q = ush.coeffQ(N);
        for(int j=0; j<N; j++) {
            System.out.println(j+" "+q[j]);
        }

        System.out.println();
        N = 10;
        double[] rho = new double[20];
        for(int i=0; i<rho.length; i++) {
            rho[i] = (i+1)*0.5/rho.length;
        }
        double[] p = ush.pCalc(N, rho, T);
        for(int j=0; j<rho.length; j++) {
            System.out.println("{"+rho[j]+","+p[j]+"},");
        }
    }
}
