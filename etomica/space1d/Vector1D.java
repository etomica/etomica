package etomica.space1d;

import etomica.space.IVector;
import etomica.space.IVectorRandom;
import etomica.util.Function;
import etomica.util.IRandom;

/**
 * Implementation of the Vector class for a 1-dimensional space. In this case the vector
 * is a trivial object formed from just one element.
 */
public final class Vector1D implements IVectorRandom, java.io.Serializable { 
    
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

    public boolean equals(IVector v) {
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

    public double x(int i) {
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

    public void Ea1Tv1(double a1, IVector u) {
        Vector1D u1 = (Vector1D) u;
        x = a1 * u1.x;
    }

    public void PEa1Tv1(double a1, IVector u) {
        x += a1 * ((Vector1D) u).x;
    }

    public void PE(double a) {
        x += a;
    }

    public void TE(double a) {
        x *= a;
    }

    public void Ev1Pv2(IVector u1, IVector u2) {
        x = ((Vector1D) u1).x + ((Vector1D) u2).x;
    }

    public void Ev1Mv2(IVector u1, IVector u2) {
        Vector1D v1 = (Vector1D) u1;
        Vector1D v2 = (Vector1D) u2;
        x = v1.x - v2.x;
    }

    public double Mv1Squared(IVector u1) {
        double dx = x - ((Vector1D) u1).x;
        return dx * dx;
    }
    
    public void mod(IVector u) {
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
        x = 1.0;
    }

    public void setRandomCube(IRandom random) {
        x = random.nextDouble() - 0.5;
    }

    public void setRandomInSphere(IRandom random) {
        x = random.nextDouble() - 0.5;
    }

    public void setRandomSphere(IRandom random) {
        randomDirection(random);
    }

    public void randomDirection(IRandom random) {
        x = random.nextInt(2) * 2 - 1;
    }

    public void E(IVector u) {
        x = ((Vector1D) u).x;
    }

    public void PE(IVector u) {
        x += ((Vector1D) u).x;
    }

    public void ME(IVector u) {
        x -= ((Vector1D) u).x;
    }

    public void TE(IVector u) {
        x *= ((Vector1D) u).x;
    }

    public void DE(IVector u) {
        x /= ((Vector1D) u).x;
    }

    public double dot(IVector u) {
        return ((Vector1D) u).x * x;
    }

    public boolean isNaN() {
        return Double.isNaN(x);
    }

    public void map(Function function) {
        x = function.f(x);
    }

}