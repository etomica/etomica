/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.math.numerical;

import etomica.util.random.IRandom;

/**
 * Optimization algorithm which attempts to smooth data without moving any
 * point far from its original location.
 * 
 * @author Andrew Schultz
 */
public class AkimaSplineSmoother {

    public AkimaSplineSmoother(IRandom random) {
        m = new double[0];
        t = new double[0];
        y = new double[0];
        trialP = new double[]{0.5, 0.75};
        this.random = random;
        nTrial = 6;
        oy = new double[1];
    }

    public boolean doRandStep() {
        int idata = random.nextInt(x.length);
        int bestidy = 0;
        double bestErr = Double.POSITIVE_INFINITY;
        for (int idy=-nTrial; idy<nTrial+1; idy++) {
            if (idy == -nTrial) {
                y[idata] -= nTrial*ss[0][idata];
            }
            else {
                y[idata] += ss[0][idata];
            }
            double sumSqErr = calcErr(idata, idata);
            if (sumSqErr < bestErr) {
                bestidy = idy;
                bestErr = sumSqErr;
            }
        }
        y[idata] += (bestidy-nTrial)*ss[0][idata];
        nAttempt[0][idata]++;
        if (bestidy != 0) {
            nAccepted[0][idata]++;
        }
        if (nAttempt[0][idata] > 10) {
            pac[0][idata] = pac[0][idata]*0.9 + ((double)nAccepted[0][idata])/nAttempt[0][idata]*0.1;
            if (pac[0][idata] > 0.1) {
                ss[0][idata] *= 1.1;
            }
            else if (pac[0][idata] < 0.05) {
                ss[0][idata] /= 1.1;
            }
            nAccepted[0][idata] = 0;
            nAttempt[0][idata] = 0;
        }
        return bestidy != 0;
    }

    public boolean doInterpolateHole() {
        int N = x.length;
        int idata = 1+random.nextInt(N-2);
        int nmov = 0;
        while (idata-nmov > 1 && idata+nmov < N-2 && N-1-nmov*2 > 5 && random.nextDouble() < 0.7) {
            nmov++;
        }
        if (oy.length < nmov*2+1) {
            oy = new double[nmov*2+1];
        }
        
        // idataMin and Max are the first and last data point that are changed
        int idataMin = idata-nmov;
        int idataMax = idata+nmov;
        for (int i=0; i<2*nmov+1; i++) {
            oy[i] = y[idataMin+i];
        }
        double oldErr = calcErr(idataMin, idataMax);

        int holeSize = idataMax-idataMin+1;
        int lN = N - holeSize;

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
        if (maxm > N-1) maxm = N-1;
        
        int j = min;
        for (int i=min; i<max+1; i++) {
            if (i<idataMin || i > idataMax) {
                lx[j] = x[i];
                ly[j] = y[i];
                j++;
            }
        }
        
        
        for (int i=minm; i<maxm-holeSize+1; i++) {
            m[i] = (ly[i+1]-ly[i])/(lx[i+1]-lx[i]);
        }
        // special-case first t
        t[0] = 0.5 * (3*m[0] - m[1]);
        for (int i=mint; i<maxt-holeSize+1; i++) {
            double am2m1;
            if (i>1) {
                // special-case second t
                am2m1 = Math.abs(m[i-1]-m[i-2]);
            }
            else {
                am2m1 = Math.abs(m[1] - m[0]);
            }
            double am4m3;
            if (i<lN-2) {
                am4m3 = Math.abs(m[i+1]-m[i]);
            }
            else {
                // special-case next-to-last t
                am4m3 = Math.abs(m[lN-2] - m[lN-3]);
            }
            if (am2m1 != am4m3) {
                t[i] = (am4m3*m[i-1] + am2m1*m[i])/(am4m3 + am2m1);
            }
            else {
                t[i] = 0.5 * (m[i-1] + m[i]);
            }
        }
        // special-case last t
        t[lN-1] = 0.5 * (3*m[lN-2]-m[lN-3]);
        
        double p0 = y[idataMin-1];
        double p1 = t[idataMin-1];
        double dx = x[idataMax+1] - x[idataMin-1];
        double p2 = (3*m[idataMin-1] - 2*t[idataMin-1] - t[idataMin])/dx;
        double p3 = (-2*m[idataMin-1] + t[idataMin-1] + t[idataMin])/(dx*dx);
        
        for (int i=0; i<idataMax-idataMin+1; i++) {
            dx = x[idataMin+i]-x[idataMin-1];
            y[idataMin+i] = p0 + (p1 + (p2 + p3*dx)*dx)*dx;
        }
        
        double newErr = calcErr(idataMin, idataMax);

        if (newErr > oldErr) {
            for (int i=0; i<2*nmov+1; i++) {
                y[idataMin+i] = oy[i];
            }
            return false;
        }        
        return true;
    }
    
    public boolean doInterpolateHole2() {
        int N = x.length;
        int idata = 2+random.nextInt(N-4);
        int nmov = 1;
        while (idata-nmov > 1 && idata+nmov < N-2 && random.nextDouble() < 0.7) {
            nmov++;
        }
        if (oy.length < nmov*2+1) {
            oy = new double[nmov*2+1];
        }
        
        // idataMin and Max are the first and last data point that are changed
        int idataMin = idata-nmov;
        int idataMax = idata+nmov;
        for (int i=0; i<2*nmov+1; i++) {
            oy[i] = y[idataMin+i];
        }
        double oldErr = calcErr(idataMin, idataMax);
        y[idata] += 0.5*(random.nextDouble()-0.5)*ss[2][idata];

        int holeSize = nmov*2;
        int lN = N - holeSize;

        int min = 0;
        int max = N-1;
        if (idataMin > 3) min = idataMin-3;
        if (idataMax < N-4) max = idataMax + 3;
        
        int j = min;
        for (int i=min; i<max+1; i++) {
            if (i<idataMin || i > idataMax || i == idata) {
                lx[j] = x[i];
                ly[j] = y[i];
                j++;
            }
        }
        
        int mint = idataMin-1;
        int maxt = mint+2;
        if (mint == 0) mint = 1;
        if (maxt == N-1) maxt = N-1;
        int minm = mint-2;
        if (minm<0) minm = 0;
        int maxm = maxt+1;
        if (maxm > lN-2) maxm = lN-2;
        
        for (int i=minm; i<maxm-holeSize+1; i++) {
            m[i] = (ly[i+1]-ly[i])/(lx[i+1]-lx[i]);
        }

        // special-case first t
        t[0] = 0.5 * (3*m[0] - m[1]);
        for (int i=mint; i<maxt+1; i++) {
            double am2m1;
            if (i>1) {
                // special-case second t
                am2m1 = Math.abs(m[i-1]-m[i-2]);
            }
            else {
                am2m1 = Math.abs(m[1] - m[0]);
            }
            double am4m3;
            if (i<lN-2) {
                am4m3 = Math.abs(m[i+1]-m[i]);
            }
            else {
                // special-case next-to-last t
                am4m3 = Math.abs(m[lN-2] - m[lN-3]);
            }
            if (am2m1 != am4m3) {
                t[i] = (am4m3*m[i-1] + am2m1*m[i])/(am4m3 + am2m1);
            }
            else {
                t[i] = 0.5 * (m[i-1] + m[i]);
            }
        }
        // special-case last t
        t[lN-1] = 0.5 * (3*m[lN-2]-m[lN-3]);

        j = idataMin-1;
        double p0 = y[j];
        double p1 = t[j];
        double dx = lx[j+1] - lx[j];
        double p2 = (3*m[j] - 2*t[j] - t[j+1])/dx;
        double p3 = (-2*m[j] + t[j] + t[j+1])/(dx*dx);
        
        for (int i=0; i<nmov; i++) {
            dx = x[idataMin+i]-x[idataMin-1];
            y[idataMin+i] = p0 + (p1 + (p2 + p3*dx)*dx)*dx;
        }

        j++;
        p0 = y[idata];
        p1 = t[idata];
        dx = lx[j+1] - lx[j];
        p2 = (3*m[j] - 2*t[j] - t[j+1])/dx;
        p3 = (-2*m[j] + t[j] + t[j+1])/(dx*dx);
        
        for (int i=0; i<nmov; i++) {
            dx = x[idata+i]-x[idata];
            y[idata+i] = p0 + (p1 + (p2 + p3*dx)*dx)*dx;
        }
        
        double newErr = calcErr(idataMin, idataMax);

        nAttempt[2][idata]++;
        boolean accepted = true;

        if (newErr > oldErr) {
            for (int i=0; i<2*nmov+1; i++) {
                y[idataMin+i] = oy[i];
            }
            accepted = false;
        }
        else {
            nAccepted[2][idata]++;
        }

        if (nAttempt[2][idata] > 10) {
            pac[2][idata] = pac[2][idata]*0.9 + ((double)nAccepted[2][idata])/nAttempt[2][idata]*0.1;
            if (pac[2][idata] > 0.1) {
                ss[2][idata] *= 1.1;
            }
            else if (pac[2][idata] < 0.05) {
                ss[2][idata] /= 1.1;
            }
            nAccepted[2][idata] = 0;
            nAttempt[2][idata] = 0;
        }

        return accepted;
    }
    
    public void doStep() {
        for (int i=0; i<100000; i++) {
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
                default:
                    throw new RuntimeException("oops");
            }
        }
        totSqErr = calcErr(0, x.length-1);
    }
    
    protected double calcErr(int idataMin, int idataMax) {
        int N = x.length;
        
        sumSqDy = 0;
        for (int idata = idataMin; idata<idataMax+1; idata++) {
            if (ey[idata] > 0) {
                double dy = (y[idata] - y0[idata]) / ey[idata];
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
                double dx2 = 0.5*(x[i+1]-x[i-1]);
                double d = (d2-lastd2);
                sumSqD2D += dx2*d*d;
                d = d2next-d2;
                sumSqD2 += dx*(d2*d2 + d2*d + d*d);
                d = d3-lastd3;
                sumSqD3D += dx2*d*d;
            }
            sumSqD3 += dx*d3*d3;
            lastd2 = d2next;
            lastd3 = d3;
        }
        return sumSqDy + d2fac*sumSqD2 + d2dfac*sumSqD2D + d3fac*sumSqD3 + d3dfac*sumSqD3D;
    }
    
    public void setSmoothedY(double[] smoothedY) {
        System.arraycopy(smoothedY, 0, y, 0, smoothedY.length);
    }

    public double[][] getDy12(int nSubPoints) {
        if (dy12 == null || dy12[0].length != (y.length-1)*nSubPoints+1) {
            dy12 = new double[2][(y.length-1)*nSubPoints+1];
        }
        int N = x.length;
        if (t.length != N) {
            m = new double[N-1];
            t = new double[N];
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
            double p1 = t[i];
            double dx = (x[i+1]-x[i]);
            double p2 = (3*m[i] - 2*t[i] - t[i+1])/dx;
            double p3 = (-2*m[i] + t[i] + t[i+1])/(dx*dx);
            for (int j=0; j<nSubPoints; j++) {
                dx = j*(x[i+1] - x[i])/nSubPoints;
                dy12[0][i*nSubPoints+j] = p1 + (2*p2 + 3*p3*dx)*dx;
                dy12[1][i*nSubPoints+j] = 2*p2 + 6*p3*dx;
            }
            if (i==N-2) {
                dx = (x[i+1] - x[i]);
                dy12[0][(N-1)*nSubPoints] = p1 + (2*p2 + 3*p3*dx)*dx;
                dy12[1][(N-1)*nSubPoints] = 2*p2 + 6*p3*dx;
            }
        }
        return dy12;
    }

    public void setInputData(double[] originalX, double[] originalY, double[] ey) {
        this.x = originalX;
        this.y0 = originalY;
        this.ey = ey;
        int N = originalX.length;
        m = new double[N-1];
        t = new double[N];
        y = new double[N];
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
                    int ii=i-1;
                    while (ey[ii] >= Double.MAX_VALUE) {
                        ii--;
                    }
                    ss[j][i] = 0.5*ey[ii];
                }
                pac[j][i] = 0.1;
            }
        }
        lx = new double[N];
        ly = new double[N];
    }

    public double getD2fac() {
        return d2fac;
    }

    /**
     * sets weight given to minimizing the 2nd derivative
     */
    public void setD2fac(double d2fac) {
        this.d2fac = d2fac;
    }

    public double getD2dfac() {
        return d2dfac;
    }

    /**
     * sets weight given to minimizing discontinuities in the 2nd derivative
     * (at data points)
     */
    public void setD2dfac(double d2dfac) {
        this.d2dfac = d2dfac;
    }

    public double getD3fac() {
        return d3fac;
    }

    /**
     * sets weight given to minimizing the 3nd derivative
     */
    public void setD3fac(double d3fac) {
        this.d3fac = d3fac;
    }

    /**
     * sets weight given to minimizing discontinuities in the 3nd derivative
     * (at data points)
     */
    public double getD3dfac() {
        return d3dfac;
    }

    public void setD3dfac(double d3dfac) {
        this.d3dfac = d3dfac;
    }

    protected double[] x, y0, ey;
    protected double[] m, t, y;
    protected double[] lx, ly;
    protected double[][] ss, pac;
    protected double[][] dy12;
    protected double totSqErr, sumSqDy, sumSqD2, sumSqD2D, sumSqD3, sumSqD3D;
    protected double d2fac, d2dfac, d3fac, d3dfac;
    protected double[] trialP;
    protected IRandom random;
    protected int nTrial;
    protected int[][] nAttempt, nAccepted;
    protected double[] oy;
}
