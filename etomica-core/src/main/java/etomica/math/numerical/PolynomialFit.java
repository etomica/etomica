/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.math.numerical;

import Jama.Matrix;

/**
 * Class that performs a polynomial fit to x,y data, optionally taking weights
 * associated with each data point.
 * 
 * @author Andrew Schultz
 */
public class PolynomialFit {

    /**
     * Perform polynomial fit of the given order.  order=2 would fit
     *     y = a*x^2 + b*x + c
     * The polynomial coefficients are returned as a double array in order of
     * increasing power of x.  So,
     *   double[] r = doFit();
     * r[i] is the coefficient for x^i
     */
    public static double[] doFit(int order, double[] x, double[] y) {
        double[] w = new double[x.length];
        for (int i=0; i<x.length; i++) {
            w[i] = 1;
        }
        return doFit(order, x, y, w);
    }

    /**
     * Perform polynomial fit of the given order with the given weights (w).
     * order=2 would fit
     *     y = a*x^2 + b*x + c
     * The polynomial coefficients are returned as a double array in order of
     * increasing power of x.  So,
     *   double[] r = doFit();
     * r[i] is the coefficient for x^i
     */
    public static double[] doFit(int order, double[] x, double[] y, double[] w) {
        if (x.length != y.length || x.length != w.length || x.length < order+1) {
            // We need at least order+1 data points to do a meaningful fit.
            return null;
        }

        double[][] M = new double[order+1][order+1];
        double[] b = new double[order+1];
        for (int i=0; i<x.length; i++) {
            if (w[i] == 0) continue;
            double xp = w[i];
            for (int ipower = 0; ipower<2*order+1; ipower++) {
                for (int irow=order; irow>-1; irow--) {
                    int col = ipower - irow;
                    if (col > -1 && col < order+1) {
                        M[irow][col] += xp;
                    }
                }
                if (ipower <= order) {
                    b[ipower] += xp * y[i];
                }
                xp *= x[i];
            }
        }
        Matrix mat = new Matrix(M);
        Matrix sol = mat.solve(new Matrix(b,order+1));
        double[] result = new double[order+1];
        for (int i=0; i<result.length; i++) {
            result[i] = sol.get(i, 0);
        }
        return result;
    }

    /**
     * Returns average relative deviation (relative to the uncertainty) of the
     * given fit from the given data.
     * <p>
     * chi = (<(diff/error)^2>/N)^.5
     */
    public static double getChi(double[] x, double[] y, double[] w, double[] poly) {
        double chiSqSum = 0, nData = 0;
        for (int i = 0; i < x.length; i++) {
            if (w[i] == 0) continue;
            double xp = 1;
            double yp = 0;
            for (int ipower = 0; ipower < poly.length; ipower++) {
                yp += poly[ipower] * xp;
                xp *= x[i];
            }
            double diff = y[i] - yp;
            chiSqSum += diff * diff * w[i];
            nData++;
        }
        return Math.sqrt(chiSqSum / nData);
    }
    
    public static void main(String[] args) {
        double[] c = doFit(2, new double[]{1,2,3}, new double[]{9,6,1});
        System.out.println("a="+c[2]+" b="+c[1]+" c="+c[0]);
    }
}
