/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.space1d;

import etomica.math.function.IFunction;
import etomica.util.random.IRandom;
import etomica.space.Vector;
import etomica.exception.MethodNotImplementedException;

/**
 * Implementation of the Vector class for a 1-dimensional space. In this case the vector
 * is a trivial object formed from just one element.
 */
public final class Vector1D implements Vector, java.io.Serializable {
    
    double x;
    private static final long serialVersionUID = 1L;

    public Vector1D() {
        x = 0.0;
    }

    public Vector1D(double a1) {
        x = a1;
    }

    public Vector1D(double[] a) {
        x = a[0];
    }

    public Vector1D(Vector1D u) {
        this.E(u);
    }

    public boolean equals(Vector v) {
        return (x == ((Vector1D) v).x);
    }

    public boolean isZero() {
        return x == 0;
    }

    public String toString() {
        return "(" + x + ")";
    }

    public int getD() {
        return 1;
    }

    public double getX(int i) {
        return x;
    }

    public void setX(int i, double d) {
        x = d;
    }

    public void assignTo(double[] array) {
        array[0] = x;
    }

    public void E(double a) {
        x = a;
    }

    public void E(double[] a) {
        x = a[0];
    }

    public void Ea1Tv1(double a1, Vector u) {
        Vector1D u1 = (Vector1D) u;
        x = a1 * u1.x;
    }

    public void PEa1Tv1(double a1, Vector u) {
        x += a1 * ((Vector1D) u).x;
    }

    public void PE(double a) {
        x += a;
    }

    public void TE(double a) {
        x *= a;
    }

    public void Ev1Pv2(Vector u1, Vector u2) {
        x = ((Vector1D) u1).x + ((Vector1D) u2).x;
    }

    public void Ev1Mv2(Vector u1, Vector u2) {
        Vector1D v1 = (Vector1D) u1;
        Vector1D v2 = (Vector1D) u2;
        x = v1.x - v2.x;
    }

    public double Mv1Squared(Vector u1) {
        double dx = x - ((Vector1D) u1).x;
        return dx * dx;
    }
    
    public void mod(Vector u) {
        mod((Vector1D) u);
    }

    public void mod(Vector1D u) {
        while (x > u.x)
            x -= u.x;
        while (x < 0.0)
            x += u.x;
    }

    public double squared() {
        return x * x;
    }

    public void normalize() {
        x = x > 0 ? 1 : -1;
    }

    public void setRandomCube(IRandom random) {
        x = random.nextFixedDouble() - 0.5;
    }

    public void setRandomInSphere(IRandom random) {
        x = 2*(random.nextFixedDouble() - 0.5);
    }

    @Override
    public void nearestImage(Vector dimensions) {
        Vector1D dimensions1D = ((Vector1D) dimensions);
        final double halfX = dimensions1D.x / 2;

        while (x > halfX)
            x -= dimensions1D.x;
        while (x < -halfX)
            x += dimensions1D.x;
    }

    public void setRandomSphere(IRandom random) {
        randomDirection(random);
    }

    public void randomDirection(IRandom random) {
        x = random.nextInt(2) * 2 - 1;
    }

    public void E(Vector u) {
        x = ((Vector1D) u).x;
    }

    public void PE(Vector u) {
        x += ((Vector1D) u).x;
    }

    public void ME(Vector u) {
        x -= ((Vector1D) u).x;
    }

    public void TE(Vector u) {
        x *= ((Vector1D) u).x;
    }

    public void DE(Vector u) {
        x /= ((Vector1D) u).x;
    }

    public double dot(Vector u) {
        return ((Vector1D) u).x * x;
    }

    public boolean isNaN() {
        return Double.isNaN(x);
    }

    public void map(IFunction function) {
        x = function.f(x);
    }

    public void XE(Vector u) {
        throw new MethodNotImplementedException();
    }
}
