package etomica.space2d;

import etomica.api.IFunction;
import etomica.api.IRandom;
import etomica.api.IVector;
import etomica.space.IVectorRandom;

/**
 * Implementation of the Vector class for a 2-dimensional space.
 */
public class Vector2D implements IVectorRandom, java.io.Serializable {

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

    public boolean equals(IVector v) {
        return (x == ((Vector2D) v).x) && (y == ((Vector2D) v).y);
    }

    public boolean isZero() {
        return (x == 0.0) && (y == 0.0);
    }

    public int getD() {
        return 2;
    }

    public double x(int i) {
        return (i == 0) ? x : y;
    }

    public void setX(int i, double d) {
        if (i == 0)
            x = d;
        else
            y = d;
    }

    public void E(IVector u) {
        x = ((Vector2D) u).x;
        y = ((Vector2D) u).y;
    }

    public void E(double[] u) {
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

    public void Ea1Tv1(double a1, IVector u) {
        Vector2D u1 = (Vector2D) u;
        x = a1 * u1.x;
        y = a1 * u1.y;
    }

    public void PEa1Tv1(double a1, IVector u) {
        Vector2D u1 = (Vector2D) u;
        x += a1 * u1.x;
        y += a1 * u1.y;
    }

    public void PE(IVector u) {
        x += ((Vector2D) u).x;
        y += ((Vector2D) u).y;
    }

    public void PE(double a) {
        x += a;
        y += a;
    }

    public void ME(IVector u) {
        x -= ((Vector2D) u).x;
        y -= ((Vector2D) u).y;
    }

    public void TE(double a) {
        x *= a;
        y *= a;
    }

    public void TE(IVector u) {
        x *= ((Vector2D) u).x;
        y *= ((Vector2D) u).y;
    }

    public void DE(IVector u) {
        x /= ((Vector2D) u).x;
        y /= ((Vector2D) u).y;
    }

    public double Mv1Squared(IVector u) {
        Vector2D u1 = (Vector2D) u;
        double dx = x - u1.x;
        double dy = y - u1.y;
        return dx * dx + dy * dy;
    }

    public void Ev1Pv2(IVector u1, IVector u2) {
        Vector2D v1 = (Vector2D) u1;
        Vector2D v2 = (Vector2D) u2;
        x = v1.x + v2.x;
        y = v1.y + v2.y;
    }

    public void Ev1Mv2(IVector u1, IVector u2) {
        Vector2D v1 = (Vector2D) u1;
        Vector2D v2 = (Vector2D) u2;
        x = v1.x - v2.x;
        y = v1.y - v2.y;
    }

    public void mod(IVector u) {
        Vector2D u2 = (Vector2D) u;
        while (x > u2.x)
            x -= u2.x;
        while (x < 0.0)
            x += u2.x;
        while (y > u2.y)
            y -= u2.y;
        while (y < 0.0)
            y += u2.y;
    }

    public double squared() {
        return x * x + y * y;
    }

    public double dot(IVector u) {
        return x * ((Vector2D) u).x + y * ((Vector2D) u).y;
    }

    public void normalize() {
        double norm = Math.sqrt(1 / (x * x + y * y));
        x *= norm;
        y *= norm;
    }

    public void setRandomCube(IRandom random) {
        x = random.nextDouble() - 0.5;
        y = random.nextDouble() - 0.5;
    }

    public void setRandomSphere(IRandom random) {
        x = Math.cos(2 * Math.PI * random.nextDouble());
        y = Math.sqrt(1.0 - x * x);
        if (random.nextDouble() < 0.5)
            y = -y;
    }

    public void setRandomInSphere(IRandom random) {
        double z1 = 0.0;
        double z2 = 0.0;
        double rsq;
        do {
            z1 = 1.0 - 2.0 * random.nextDouble();
            z2 = 1.0 - 2.0 * random.nextDouble();
            rsq = z1 * z1 + z2 * z2;
        } while (rsq > 1.0);
        x = z1;
        y = z2;
    }

    public boolean isNaN() {
        return Double.isNaN(x) || Double.isNaN(y);
    }

    public void map(IFunction function) {
        x = function.f(x);
        y = function.f(y);
    }

}