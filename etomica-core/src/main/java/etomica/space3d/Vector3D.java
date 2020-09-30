/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.space3d;

import etomica.math.function.IFunction;
import etomica.space.Vector;
import etomica.util.random.IRandom;

import java.util.Objects;

/**
 * Implementation of the Vector class for a 3-dimensional space.
 */
public final class Vector3D implements Vector, java.io.Serializable {

    protected double x, y, z;
    private static final long serialVersionUID = 1L;

    public int getD() {
        return 3;
    }

    public Vector3D() {
        x = 0.0;
        y = 0.0;
        z = 0.0;
    }

    public Vector3D(double a1, double a2, double a3) {
        x = a1;
        y = a2;
        z = a3;
    }

    public Vector3D(double[] a) {
        if(a.length != 3){ 
            throw new IllegalArgumentException("Vector3D must be given a 3 element array.");
        }
        x = a[0];
        y = a[1];
        z = a[2];
    }

    public Vector3D(Vector3D u) {
        this.E(u);
    }

    public String toString() {
        return "(" + x + ", " + y + ", " + z + ")";
    }

    public double getX(int i) {
        return ((i == 0) ? x : (i == 1) ? y : z);
    }

    public void assignTo(double[] array) {
        array[0] = x;
        array[1] = y;
        array[2] = z;
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null) return false;
        if (!(o instanceof Vector) || ((Vector) o).getD() != this.getD()) {
            return false;
        }
        Vector v = ((Vector) o);
        return v.x() == this.x() && v.y() == this.y() && v.z() == this.z();
    }

    public boolean equals(Vector v) {
        return (x == v.x()) && (y == v.y())
                && (z == v.z());
    }

    public boolean isZero() {
        return (x == 0.0) && (y == 0.0) && (z == 0);
    }

    public void E(Vector u) {
        x = u.getX(0);
        y = u.getX(1);
        z = u.getX(2);
    }

    public void E(double a) {
        x = a;
        y = a;
        z = a;
    }

    public void E(double a, double b, double c) {
        x = a;
        y = b;
        z = c;
    }

    public void E(double... u) {
        if(u.length != 3){ 
            throw new IllegalArgumentException("Vector3D must be given a 3 element array.");
        }
        x = u[0];
        y = u[1];
        z = u[2];
    }
    
    public void Ea1Tv1(double a1, Vector u) {
        x = a1 * u.getX(0);
        y = a1 * u.getX(1);
        z = a1 * u.getX(2);
    }

    public void PEa1Tv1(double a1, Vector u) {
        x += a1 * u.getX(0);
        y += a1 * u.getX(1);
        z += a1 * u.getX(2);
    }

    public void PE(Vector u) {
        x += u.getX(0);
        y += u.getX(1);
        z += u.getX(2);
    }

    public void PE(double a) {
        x += a;
        y += a;
        z += a;
    }

    public void ME(Vector u) {
        x -= u.getX(0);
        y -= u.getX(1);
        z -= u.getX(2);
    }

    public void TE(double a) {
        x *= a;
        y *= a;
        z *= a;
    }

    public void TE(Vector u) {
        x *= u.getX(0);
        y *= u.getX(1);
        z *= u.getX(2);
    }

    public void DE(Vector u) {
        x /= u.getX(0);
        y /= u.getX(1);
        z /= u.getX(2);
    }

    public void Ev1Pv2(Vector v1, Vector v2) {
        x = v1.getX(0) + v2.getX(0);
        y = v1.getX(1) + v2.getX(1);
        z = v1.getX(2) + v2.getX(2);
    }

    public void Ev1Mv2(Vector v1, Vector v2) {
        x = v1.getX(0) - v2.getX(0);
        y = v1.getX(1) - v2.getX(1);
        z = v1.getX(2) - v2.getX(2);
    }

    public void mod(Vector v) {
        while (x > v.x())
            x -= v.x();
        while (x < 0.0)
            x += v.x();
        while (y > v.y())
            y -= v.y();
        while (y < 0.0)
            y += v.y();
        while (z > v.z())
            z -= v.z();
        while (z < 0.0)
            z += v.z();
    }

    public double squared() {
        return x * x + y * y + z * z;
    }

    public double Mv1Squared(Vector v1) {
        double dx = x - v1.getX(0);
        double dy = y - v1.getX(1);
        double dz = z - v1.getX(2);
        return dx * dx + dy * dy + dz * dz;
    }

    public double dot(Vector u) {
        return x * u.getX(0) + y * u.getX(1) + z * u.getX(2);
    }

    /*
     * Sets this vector to an arbitrary vector in the plane normal to the given vector.
     * Does not normalize this vector on completion.  Assumes that the given vector is not identically zero.
     */
    public void setPerpendicularTo(Vector v) {
        Vector3D v3 = (Vector3D)v;
        if(v3.z != 0) {
            x = v3.z;
            y = v3.z;
            z = -(v3.x + v3.y);
        } else {
            x = -(v3.y + v3.z);
            y = v3.x;
            z = v3.x;
        }
    }


    public void setX(int a, double d) {
        if (a == 0)
            x = d;
        else if (a == 1)
            y = d;
        else
            z = d;
    }

    public void setRandomCube(IRandom random) {
        x = random.nextFixedDouble() - 0.5;
        y = random.nextFixedDouble() - 0.5;
        z = random.nextFixedDouble() - 0.5;
    }

    //generate point on surface of sphere according to method (4) here:
    //http://www.math.niu.edu/~rusin/known-math/96/sph.rand
    //and scale into the interior according to r^2
    public void setRandomInSphere(IRandom random) {
        double r = cubeRoot(random.nextFixedDouble());
        double u, v, s;
        do {
            u = 1.0 - 2.0*random.nextFixedDouble();
            v = 1.0 - 2.0*random.nextFixedDouble();
            s = u*u + v*v;
        } while(s > 1);
        double ra = 2.*r * Math.sqrt(1.-s);
        x = ra * u;
        y = ra * v;
        z = r * (2*s- 1.);
    }

    @Override
    public void nearestImage(Vector dim) {
        final double halfX = dim.x() / 2;
        final double halfY = dim.y() / 2;
        final double halfZ = dim.z() / 2;

        while (x > halfX)
            x -= dim.x();
        while (x < -halfX)
            x += dim.x();
        while (y > halfY)
            y -= dim.y();
        while (y < -halfY)
            y += dim.y();
        while (z > halfZ)
            z -= dim.z();
        while (z < -halfZ)
            z += dim.z();
    }

    /**
     * Creating a random unit vector on unit sphere Uses only two random number
     * generator at a time
     * <p>
     * Based on M.P. Allen and D.J. Tildesley, Computer Simulation of Liquids, p 349.
     * (taken from Marsaglia, G. (1972), Choosing a point from the surface of a sphere,
     * Ann. math. Stat. 43, 645-646.)
     */
    public void setRandomSphere(IRandom random) {
        double z1, z2, zsq;
        do  {
            z1 = 2.0 * random.nextFixedDouble() - 1.0;
            z2 = 2.0 * random.nextFixedDouble() - 1.0;
            zsq = z1 * z1 + z2 * z2;
        } while (zsq > 1.0);

        double ranh = 2.0 * Math.sqrt(1.0 - zsq);
        x = z1 * ranh;
        y = z2 * ranh;
        z = 1.0 - 2.0 * zsq;
    }

    public void XE(Vector u) {//cross product
        double xNew = y * u.getX(2) - z * u.getX(1);
        double yNew = z * u.getX(0) - x * u.getX(2);
        z = x * u.getX(1) - y * u.getX(0);
        y = yNew;
        x = xNew;
    }

    public void normalize() {
        double norm = 1. / Math.sqrt(x * x + y * y + z * z);
        x *= norm;
        y *= norm;
        z *= norm;
    }

    public boolean isNaN() {
        return Double.isNaN(x) || Double.isNaN(y) || Double.isNaN(z);
    }

    public void map(IFunction function) {
        x = function.f(x);
        y = function.f(y);
        z = function.f(z);
    }

    protected static final int nLookUp = 5000;
    protected static final double[] lookUp;
    static {
        lookUp = new double[nLookUp];
        for(int i=0; i<nLookUp; i++) {
            lookUp[i] = Math.pow((double)(i+1)/nLookUp,1.0/3.0);
        }
    }

    //cube root via look-up and Newton/Halley iteration
    //see here: http://metamerist.com/cbrt/cbrt.htm
    public static double cubeRoot(double R) {
        int iR = (int) (R * nLookUp);
        // method is accurate to machine precision except for very small r.
        // fall back to pow
        if (iR < 10) return Math.pow(R, 1.0 / 3.0);
        double a = lookUp[iR];
        for (int i = 0; i < 2; i++) {
            double a3 = a * a * a;
            double t1 = a3 + R;
            a = a * (t1 + R) / (t1 + a3);
        }
        return a;
    }

    public Vector3D duplicate() {
        return new Vector3D(this.x, this.y, this.z);
    }

    @Override
    public int hashCode() {
        return Objects.hash(x, y, z);
    }
}
