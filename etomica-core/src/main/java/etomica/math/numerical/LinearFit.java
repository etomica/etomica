/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.math.numerical;

import Jama.Matrix;

/**
 * Class that performs a linear fit to x,y data, optionally taking weights
 * associated with each data point.
 * 
 * @author Andrew Schultz
 */
public class LinearFit {

    public static double[] doFit(double[] x, double[] y) {
        double[] w = new double[x.length];
        for (int i=0; i<x.length; i++) {
            w[i] = 1;
        }
        return doFit(x, y, w);
    }
    
    public static double[] doFit(double[] x, double[] y, double[] w) {
        if (x.length != y.length || x.length != w.length || x.length < 2) {
            // We need at least two data points to do a meaningful fit.
            return null;
        }
        double[][] M = new double[2][2];
        double[] b = new double[2];
        for (int i=0; i<x.length; i++) {
            M[0][0] += x[i]*x[i]*w[i];
            M[1][0] += x[i]*w[i];
            M[0][1] += x[i]*w[i];
            b[0] += y[i]*x[i]*w[i];
            b[1] += y[i]*w[i];
            M[1][1] += w[i];
        }
        Matrix mat = new Matrix(M);
        Matrix sol = mat.solve(new Matrix(b,2));
        double[] result = new double[2];
        result[1] = sol.get(0,0);
        result[0] = sol.get(1,0);
        return result;
    }

    public static void main(String[] args) {
        double[] r = doFit(new double[]{1,2,3}, new double[]{0,1,2});
        System.out.println("m="+r[1]+" b="+r[0]);
    }
}
