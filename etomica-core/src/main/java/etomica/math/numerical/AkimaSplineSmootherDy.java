/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.math.numerical;

import etomica.util.random.IRandom;

/**
 * Optimization algorithm which attempts to smooth data without moving any
 * point far from its original location.  AkimaSplineSmootherDy specifically
 * attempts to find a smooth function for dy/dx, which when integrated yields a
 * function that matches the original y(x) within the original error bars.
 * 
 * @author Andrew Schultz
 */
public class AkimaSplineSmootherDy extends AkimaSplineSmoother {

    public AkimaSplineSmootherDy(IRandom random) {
        super(random);
        iy = new double[0];
        trialP = new double[]{0.5, 0.75, 0.99};
        ssInit = 0.01;
    }

    public void doShiftInitY() {
        // attempt to change the value of the smoothed original function
        // this value (iy[0]) is used as the starting point for integration
        double oldInitY = iy[0];
        double oldErr = calcErr(0, x.length-1);
        iy[0] += 2.0*ssInit*(random.nextDouble()-0.5);
        double newErr = calcErr(0, x.length-1);
        if (oldErr < newErr) {
            iy[0] = oldInitY;
        }
    }
    
    public void doStep() {
        for (int i=0; i<10000; i++) {
            double r = random.nextDouble();
            int iSelect = -1;
            for (int j=0; j<trialP.length; j++) {
                if (r < trialP[j]) {
                    iSelect = j;
                    break;
                }
            }
            if (iSelect == -1) {
                iSelect = trialP.length;
            }
            switch(iSelect) {
                case 0: 
                    doRandStep();
                    break;
                case 1:
                    doInterpolateHole();
                    break;
                case 2:
                    doInterpolateHole2();
                    break;
                case 3:
                    doShiftInitY();
                    break;
                default:
                    throw new RuntimeException("oops");
            }
        }
        totSqErr = calcErr(0, x.length-1);
        System.out.println("hi "+totSqErr);
    }
    
    protected double calcErr(int idataMin, int idataMax) {
        // we need to integrate the whole function
        // we could actually integrate from idataMin-something
        // to N if we were really careful
        idataMin = 0;
        idataMax = x.length-1;
        int N = x.length;

        // perform the actual integration
        getIy();
        
        sumSqDy = 0;
        for (int idata = idataMin; idata<idataMax+1; idata++) {
            if (ey[idata] > 0) {
                double dy = (iy[idata] - y0[idata]) / ey[idata];
                sumSqDy += dy*dy;
            }
        }

        int min = 0;
        int max = N-1;
        if (idataMin > 4) min = idataMin-4;
        if (idataMax < N-5) max = idataMax + 4;
        
        int mint = min;
        if (mint==0) mint=1;
        int maxt = max;
        if (maxt == N-1) maxt = N-2;
        int minm = min-2;
        if (minm<0) minm = 0;
        int maxm = max+1;
        if (maxm > N-2) maxm = N-2;
        
        
        for (int i=minm; i<maxm+1; i++) {
            m[i] = (y[i+1]-y[i])/(x[i+1]-x[i]);
        }
        
        // special-case first t -- if we don't need it, oh well
        t[0] = 0.5 * (3*m[0] - m[1]);
        for (int i=mint; i<maxt+1; i++) {
            double am2m1;
            if (i>1) {
                am2m1 = Math.abs(m[i-1]-m[i-2]);
            }
            else {
                // special-case second t
                am2m1 = Math.abs(m[1] - m[0]);
            }
            double am4m3;
            if (i<N-2) {
                am4m3 = Math.abs(m[i+1]-m[i]);
            }
            else {
                // special-case next-to-last t
                am4m3 = Math.abs(m[N-2] - m[N-3]);
            }
            if (am2m1 != am4m3) {
                t[i] = (am4m3*m[i-1] + am2m1*m[i])/(am4m3 + am2m1);
            }
            else {
                t[i] = 0.5 * (m[i-1] + m[i]);
            }
        }
        // special-case last t -- if we don't need it, oh well
        t[N-1] = 0.5 * (3*m[N-2]-m[N-3]);

        sumSqD2D = 0;
        sumSqD2 = 0;
        sumSqD3 = 0;
        sumSqD3D = 0;
        double lastd3 = 0, lastd2 = 0;
        for (int i=min; i<max; i++) {
            double dx = (x[i+1]-x[i]);
            double p2 = (3*m[i] - 2*t[i] - t[i+1])/dx;
            double p3 = (-2*m[i] + t[i] + t[i+1])/(dx*dx);
            double d2 = 2*p2;
            double d2next = 2*p2 + 6*p3*(x[i+1]-x[i]);
            double d3 = 6*p3;
            if (i>min) {
                double dx2 = Math.abs(0.5*(x[i+1]-x[i-1]));
                double d = (d2-lastd2);
                sumSqD2D += dx2*d*d;
                d = d2next-d2;
                sumSqD2 += dx*(d2*d2 + d2*d + d*d);
                d = d3-lastd3;
                sumSqD3D += dx2*d*d;
            }
            sumSqD3 += Math.abs(dx)*d3*d3;
            lastd2 = d2next;
            lastd3 = d3;
        }
        return sumSqDy + d2fac*sumSqD2 + d2dfac*sumSqD2D + d3fac*sumSqD3 + d3dfac*sumSqD3D;
    }

    /**
     * calculate the cumulative integral of y, interpolating dy/dx with a
     * spline curve
     */
    public double[] getIy() {
        int N = x.length;
        if (t.length != N) {
            m = new double[N-1];
            t = new double[N];
            iy = new double[N];
        }
        
        for (int i=0; i<N-1; i++) {
            m[i] = (y[i+1]-y[i])/(x[i+1]-x[i]);
        }

        // special-case first t
        t[0] = 0.5 * (3*m[0] - m[1]);
        for (int i=1; i<N-1; i++) {
            double am2m1;
            if (i>1) {
                // special-case second t
                am2m1 = Math.abs(m[i-1]-m[i-2]);
            }
            else {
                am2m1 = Math.abs(m[1] - m[0]);
            }
            double am4m3;
            if (i<N-2) {
                am4m3 = Math.abs(m[i+1]-m[i]);
            }
            else {
                // special-case next-to-last t
                am4m3 = Math.abs(m[N-2] - m[N-3]);
            }
            if (am2m1 != am4m3) {
                t[i] = (am4m3*m[i-1] + am2m1*m[i])/(am4m3 + am2m1);
            }
            else {
                t[i] = 0.5 * (m[i-1] + m[i]);
            }
        }
        // special-case last t
        t[N-1] = 0.5 * (3*m[N-2]-m[N-3]);

        for (int i=0; i<N-1; i++) {
            double p0 = y[i];
            double p1 = t[i];
            double dx = (x[i+1]-x[i]);
            double dx2 = dx*dx;
            double p2 = (3*m[i] - 2*t[i] - t[i+1])/dx;
            double p3 = (-2*m[i] + t[i] + t[i+1])/dx2;
            
            iy[i+1] = iy[i] + p0*dx + 0.5*p1*dx2 + p2*dx2*dx/3.0 + 0.25*p3*dx2*dx2;
        }
        return iy;
    }

    public void setInputData(double[] originalX, double[] originalY, double[] ey) {
        this.x = originalX;
        this.y0 = originalY;
        this.ey = ey;
        int N = originalX.length;
        m = new double[N-1];
        t = new double[N];
        y = new double[N];
        // pretend that y is actually y.  later, we'll store dy in y
        System.arraycopy(originalY, 0, y, 0, N);
        ss = new double[3][N];
        pac = new double[3][N];
        nAttempt = new int[3][N];
        nAccepted = new int[3][N];
        for (int j=0; j<3; j++) {
            for (int i=0; i<N; i++) {
                if (ey[i] < Double.MAX_VALUE) {
                    ss[j][i] = 0.5*ey[i];
                }
                else {
                    ss[j][i] = 0.5*ey[i-1];
                }
                pac[j][i] = 0.1;
            }
        }
        lx = new double[N];
        ly = new double[N];

        // calculate dy/dx, replace y with dy/dx
        double[] dy1 = getDy12(1)[0];
        System.arraycopy(dy1, 0, y, 0, N);

        // calculate iy now
        iy = new double[N];
        iy[0] = originalY[0];
        getIy();
    }

    /**
     * Set the nominal values for x and y and uncertainty in y (ey).  Also give
     * an initial estimate for dy/dx (initDy)
     */
    public void setInputData(double[] originalX, double[] originalY, double[] initDy, double[] ey) {
        this.x = originalX;
        this.y0 = originalY;
        this.ey = ey;
        int N = originalX.length;
        m = new double[N-1];
        t = new double[N];
        y = new double[N];
        ss = new double[3][N];
        pac = new double[3][N];
        nAttempt = new int[3][N];
        nAccepted = new int[3][N];
        for (int j=0; j<3; j++) {
            for (int i=0; i<N; i++) {
                if (ey[i] < Double.MAX_VALUE) {
                    ss[j][i] = 0.5*ey[i];
                }
                else {
                    ss[j][i] = 0.5*ey[i-1];
                }
                pac[j][i] = 0.1;
            }
        }
        lx = new double[N];
        ly = new double[N];

        System.arraycopy(initDy, 0, y, 0, N);

        // calculate iy now
        iy = new double[N];
        iy[0] = originalY[0];
        getIy();
    }

    protected double[] iy;
    protected double ssInit;
}
