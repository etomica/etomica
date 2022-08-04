/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.math.numerical;

import etomica.math.function.FunctionMultiDimensionalDifferentiable;

public class SteepestDescent {

    public SteepestDescent(FunctionMultiDimensionalDifferentiable f) {
        this.f = f;
    }
    
    public double[] minimize(double[] xGuess, double[] xStep, double tol, int maxIter) {
        return minimize(xGuess, xStep, tol, maxIter, false);
    }

    public double[] minimize(double[] xGuess, double[] xStep, double tol, int maxIter, boolean verbose) {
        final int n = xGuess.length;
        double[] der = new double[n];
        double[] dir = new double[n];
        double[] x = xGuess.clone();
        double stepSize = 1;
        double dt = 0;
        double lastVal = Double.MAX_VALUE;
        double[] lastDir = new double[x.length];
        double totalStep = 0, lastTotalStep = 0;
        int lastAvgStep = -1;
        for (int iter=0; iter<maxIter; iter++) {
            dt = 0;
            double dirdot = 0;
            for (int i=0; i<n; i++) {
                dirdot = lastDir[i]*dir[i];
            }
            if (dirdot < -0.2 && iter>5 && iter>lastAvgStep+3 && Math.random() < 0.5) {
                System.arraycopy(f.gradf(x), 0, der, 0, n);
                for (int i=0; i<n; i++) {
                    dir[i] = totalStep*dir[i] + lastTotalStep*lastDir[i];
                    dt += dir[i]*dir[i];
                }
                lastAvgStep = iter;
            }
            else {
                lastDir = dir.clone();
                System.arraycopy(f.gradf(x), 0, der, 0, n);
                for (int i=0; i<n; i++) {
                    dir[i] = -der[i] * xStep[i];
                    dt += dir[i]*dir[i];
                }
            }
            if (dt == 0) return x;
            dt = Math.sqrt(dt);
            double totalD = 0;
            for (int i=0; i<n; i++) {
                dir[i] /= dt;
                totalD += dir[i]*der[i];
            }
            double val = f.f(x);
            if (verbose) {
                if (iter > 0)
                    System.out.println(String.format("%4d   % 10.4e   % 10.4e   % 10.4e   % 10.4e  %10.4e % 7.3f", iter, val, val - lastVal, totalD, totalStep, stepSize, dirdot));
                else
                    System.out.println(String.format("%4d   % 10.4e                 % 10.4e   % 10.4e", iter, val, totalD, totalStep));
            }
            if (lastAvgStep == iter) System.out.println("averaging last directions");
//            System.out.println("dir: "+Arrays.toString(dir));
            if (lastVal - val < tol) return x;
//            System.out.println(0+" "+val+" "+totalD);
            lastVal = val;
            double newTotalD = 0;
            lastTotalStep = totalStep;
            totalStep = 0;

            while (true) {
                for (int i=0; i<n; i++) {
                    x[i] += dir[i]*stepSize;
                }
                totalStep += stepSize;
//                System.out.println("x  "+Arrays.toString(x));
                System.arraycopy(f.gradf(x), 0, der, 0, n);
                for (int i=0; i<n; i++) {
                    newTotalD += dir[i]*der[i];
                }
                val = f.f(x);
//                System.out.println(String.format("%4d  %3d    %10.4e   %10.4e     %10.4e", iter, 1, val, newTotalD, stepSize));
                
                if (newTotalD > 0 || val > lastVal) {
                    // way too far
//                    System.out.println("way too far, iter="+iter+" initial step");
//                    System.out.println("way too far, iter="+iter+" initial step "+stepSize+" "+(val-lastVal)+" "+newTotalD);
                    for (int i=0; i<n; i++) {
                        x[i] -= dir[i]*stepSize;
                    }
                    totalStep = 0;
                    
                    stepSize *= 0.1;
                    if (stepSize == 0) return x;
                    newTotalD = 0;
                    continue;
                }
                if ((val < lastVal && val > lastVal - tol) || Math.abs(stepSize) < 1e-14) return x;
                break;
            }
            double oldStepSize = stepSize;
            for (int j=0; j<4; j++) {
                double innerLastVal = val;
                double newStepSize = -newTotalD * oldStepSize / (newTotalD - totalD);
                totalD = newTotalD;
                for (int i=0; i<n; i++) {
                    x[i] += dir[i]*newStepSize;
                }
                val = f.f(x);
                while (val > innerLastVal) {
                    for (int i=0; i<n; i++) {
                        x[i] -= dir[i]*newStepSize;
                    }
                    newStepSize *= 0.1;
                    for (int i=0; i<n; i++) {
                        x[i] += dir[i]*newStepSize;
                    }
                    val = f.f(x);
//                    System.out.println((totalStep+newStepSize)+" "+val);
                    if (Math.abs(newStepSize) < 1e-14) break;
                }
                totalStep += newStepSize;
    
//                System.out.println("x  "+Arrays.toString(x));
                newTotalD = 0;
                System.arraycopy(f.gradf(x), 0, der, 0, n);
                for (int i=0; i<n; i++) {
                    newTotalD += dir[i]*der[i];
                }
//                System.out.println(totalStep+" "+val+" "+newTotalD);
                oldStepSize = newStepSize;
//                System.out.println(String.format("%4d  %3d    %10.4e   %10.4e     %10.4e", iter, j+2, val, newTotalD, newStepSize));
                if (val > innerLastVal - tol || Math.abs(newStepSize) < 1e-14) break;
                
                if (totalD < 0 && newTotalD > -0.5*totalD) {
                    System.out.println("way too far, outer="+iter+" inner="+j);
                }
            }

        }
        
        return x;
    }

    protected FunctionMultiDimensionalDifferentiable f;
}
