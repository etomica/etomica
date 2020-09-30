/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.space2d;

import etomica.exception.MethodNotImplementedException;
import etomica.math.function.IFunction;
import etomica.space.Vector;
import etomica.util.random.IRandom;

/**
 * Implementation of the Vector class for a 2-dimensional space.
 */
public final class Vector2D implements Vector, java.io.Serializable {

    double x, y;
    private static final long serialVersionUID = 1L;

    public Vector2D() {
        x = 0.0;
        y = 0.0;
    }

    public Vector2D(double x, double y) {
        this.x = x;
        this.y = y;
    }

    public Vector2D(double[] a) {
        x = a[0];
        y = a[1];
    }//should check length of a for exception

    public Vector2D(Vector2D u) {
        this.E(u);
    }

    public String toString() {
        return "(" + x + ", " + y + ")";
    }

    public void assignTo(double[] array) {
        array[0] = x;
        array[1] = y;
    }

    public boolean equals(Vector v) {
        return (x == v.x()) && (y == v.y());
    }

    public boolean isZero() {
        return (x == 0.0) && (y == 0.0);
    }

    public int getD() {
        return 2;
    }

    public double getX(int i) {
        return (i == 0) ? x : y;
    }

    public void setX(int i, double d) {
        if (i == 0)
            x = d;
        else
            y = d;
    }

    public void E(Vector u) {
        x = u.x();
        y = u.y();
    }

    public void E(double... u) {
        if(u.length != 2){throw new IllegalArgumentException("Vector2D must be given a 2 element array.");}
        x = u[0];
        y = u[1];
    }

    public void E(double a) {
        x = a;
        y = a;
    }

    public void E(double a, double b) {
        x = a;
        y = b;
    }

    public void Ea1Tv1(double a1, Vector v) {
        x = a1 * v.x();
        y = a1 * v.y();
    }

    public void PEa1Tv1(double a1, Vector v) {
        x += a1 * v.x();
        y += a1 * v.y();
    }

    public void PE(Vector v) {
        x += v.x();
        y += v.y();
    }

    public void PE(double a) {
        x += a;
        y += a;
    }

    public void ME(Vector v) {
        x -= v.x();
        y -= v.y();
    }

    public void TE(double a) {
        x *= a;
        y *= a;
    }

    public void TE(Vector v) {
        x *= v.x();
        y *= v.y();
    }

    public void DE(Vector v) {
        x /= v.x();
        y /= v.y();
    }

    public double Mv1Squared(Vector v1) {
        double dx = x - v1.x();
        double dy = y - v1.y();
        return dx * dx + dy * dy;
    }

    public void Ev1Pv2(Vector v1, Vector v2) {
        x = v1.x() + v2.x();
        y = v1.y() + v2.y();
    }

    public void Ev1Mv2(Vector v1, Vector v2) {
        x = v1.x() - v2.x();
        y = v1.y() - v2.y();
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
    }

    public double squared() {
        return x * x + y * y;
    }

    public double dot(Vector v) {
        return x * v.x() + y * v.y();
    }

    public void normalize() {
        double norm = Math.sqrt(1 / (x * x + y * y));
        x *= norm;
        y *= norm;
    }

    public void setRandomCube(IRandom random) {
        x = random.nextFixedDouble() - 0.5;
        y = random.nextFixedDouble() - 0.5;
    }

    public void setRandomSphere(IRandom random) {
        x = Math.cos(2 * Math.PI * random.nextFixedDouble());
        y = Math.sqrt(1.0 - x * x);
        if (random.nextInt(2) == 0)
            y = -y;
    }

    public void setRandomInSphere(IRandom random) {
        double z1 = 0.0;
        double z2 = 0.0;
        double rsq;
        do {
            z1 = 1.0 - 2.0 * random.nextFixedDouble();
            z2 = 1.0 - 2.0 * random.nextFixedDouble();
            rsq = z1 * z1 + z2 * z2;
        } while (rsq > 1.0);
        x = z1;
        y = z2;
    }

    @Override
    public void nearestImage(Vector dimensions) {
        final double halfX = dimensions.x() / 2;
        final double halfY = dimensions.y() / 2;

        while (x > halfX)
            x -= dimensions.x();
        while (x < -halfX)
            x += dimensions.x();
        while (y > halfY)
            y -= dimensions.y();
        while (y < -halfY)
            y += dimensions.y();
    }

    @Override
    public Vector duplicate() {
        return new Vector2D(x, y);
    }

    public boolean isNaN() {
        return Double.isNaN(x) || Double.isNaN(y);
    }

    public void map(IFunction function) {
        x = function.f(x);
        y = function.f(y);
    }

    public void XE(Vector u) {
        throw new MethodNotImplementedException();
    }
}
