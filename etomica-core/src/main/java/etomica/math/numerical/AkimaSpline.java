/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.math.numerical;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;

/**
 * DataProcessor that interpolates incoming data using an Akima spline
 * interpolation scheme.
 * 
 *   "A New Method of Interpolation and Smooth Curve Fitting Based on Local
 *   Procedures", H. Akima, JACM 17 (1970), 589
 * 
 * @author Andrew Schultz
 */
public class AkimaSpline {

    public AkimaSpline() {
        m = new double[0];
        t = new double[0];
        x = new double[0];
        y = new double[0];
        interpolatedY = new double[0];
        interpolatedDy = new double[0];
    }

    public double[] doInterpolation(double[] interpolatedX) {
        int N = x.length;
        if (t.length != N) {
            m = new double[N-1];
            t = new double[N];
        }
        if (interpolatedY.length != interpolatedX.length) {
            interpolatedY = new double[interpolatedX.length];
        }

        if (dirty) {
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
        }
        
        dirty = false;

        int j = 0;
        for (int i=0; j<interpolatedX.length && i<N-1; i++) {
            if (i < N-2 && interpolatedX[j] > x[i+1]) continue;
            double p0 = y[i];
            double p1 = t[i];
            double dx = (x[i+1]-x[i]);
            double p2 = (3*m[i] - 2*t[i] - t[i+1])/dx;
            double p3 = (-2*m[i] + t[i] + t[i+1])/(dx*dx);

            while (j < interpolatedX.length && (i == N-2 || interpolatedX[j] <= x[i+1])) {
                dx = interpolatedX[j] - x[i];
                interpolatedY[j] = p0 + p1*dx + p2*dx*dx + p3*dx*dx*dx;
                j++;
            }
        }
        return interpolatedY;
    }

    public double[] doInterpolationDy(double[] interpolatedX) {
        int N = x.length;
        if (t.length != N) {
            m = new double[N-1];
            t = new double[N];
        }
        if (interpolatedDy == null || interpolatedDy.length != interpolatedX.length) {
            interpolatedDy = new double[interpolatedX.length];
        }

        if (dirty) {
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
        }

        int j = 0;
        for (int i=0; j<interpolatedX.length && i<N-1; i++) {
            if (i < N-2 && interpolatedX[j] > x[i+1]) continue;
            double p1 = t[i];
            double dx = (x[i+1]-x[i]);
            double p2 = (3*m[i] - 2*t[i] - t[i+1])/dx;
            double p3 = (-2*m[i] + t[i] + t[i+1])/(dx*dx);

            while (j < interpolatedX.length && (i == N-2 || interpolatedX[j] <= x[i+1])) {
                dx = interpolatedX[j] - x[i];
                interpolatedDy[j] = p1 + 2*p2*dx + 3*p3*dx*dx;
                j++;
            }
        }
        return interpolatedDy;
    }

    public static void main(String[] args) {
        FileReader fileReader;
        String infile = "B7.dat.scaled";
        if (args.length > 0) {
            infile = args[0];
        }
        try {
            fileReader = new FileReader(infile);
        }catch(IOException e) {
            throw new RuntimeException("Cannot open infile, caught IOException: " + e.getMessage());
        }
        ArrayList<Double> xList = new ArrayList<Double>();
        ArrayList<Double> yList = new ArrayList<Double>();
        try {
            BufferedReader bufReader = new BufferedReader(fileReader);
            while (true) {
                String line = bufReader.readLine();
                if (line == null) {
                    break;
                }
                String[] xy = line.split(" +");
                xList.add(Double.parseDouble(xy[0]));
                yList.add(Double.parseDouble(xy[1]));
            }
            fileReader.close();
        } catch(IOException e) {
            throw new RuntimeException("Problem reading d.dat, caught IOException: " + e.getMessage());
        }

        double[] x = new double[xList.size()];
        double[] y = new double[yList.size()];
        for (int i=0; i<x.length; i++) {
            x[i] = xList.get(i);
            y[i] = yList.get(i);
        }
        AkimaSpline fitter = new AkimaSpline();
        fitter.setInputData(x, y);
        double[] ix = new double[10*(x.length-1)+1];
        for (int i=0; i<x.length-1; i++) {
            for (int j=0; j<10; j++) {
                ix[i*10+j] = x[i] + j*(x[i+1]-x[i])/10.0;
            }
        }
        ix[ix.length-1] = x[x.length-1];
        double[] iy = fitter.doInterpolation(ix);
        for (int i=0; i<ix.length; i++) {
            System.out.println(ix[i]+" "+iy[i]);
        }
    }

    public void markDirty() {
        dirty = true;
    }

    public void setInputData(double[] x, double[] y) {
        this.x = x;
        this.y = y;
        dirty = true;
    }

    protected double[] m, t, x, y, interpolatedY, interpolatedDy;
    protected boolean dirty = false;
}
