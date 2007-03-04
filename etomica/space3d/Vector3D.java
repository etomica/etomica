package etomica.space3d;

import etomica.simulation.Simulation;
import etomica.space.IVector;
import etomica.space.IVectorRandom;
import etomica.util.Function;
import etomica.util.IRandom;

/**
 * Implementation of the Vector class for a 3-dimensional space.
 */
public final class Vector3D implements IVectorRandom, java.io.Serializable {

    double x, y, z;
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

    public double x(int i) {
        return ((i == 0) ? x : (i == 1) ? y : z);
    }

    public void assignTo(double[] array) {
        array[0] = x;
        array[1] = y;
        array[2] = z;
    }

    public boolean equals(IVector v) {
        return (x == ((Vector3D) v).x) && (y == ((Vector3D) v).y)
                && (z == ((Vector3D) v).z);
    }

    public boolean isZero() {
        return (x == 0.0) && (y == 0.0) && (z == 0);
    }

    public void E(IVector u) {
        x = ((Vector3D) u).x;
        y = ((Vector3D) u).y;
        z = ((Vector3D) u).z;
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

    public void E(double[] u) {
        if(u.length != 3){ 
            throw new IllegalArgumentException("Vector3D must be given a 3 element array.");
        }
        x = u[0];
        y = u[1];
        z = u[2];
    }
    
    public void Ea1Tv1(double a1, IVector u) {
        x = a1 * ((Vector3D) u).x;
        y = a1 * ((Vector3D) u).y;
        z = a1 * ((Vector3D) u).z;
    }

    public void PEa1Tv1(double a1, IVector u) {
        x += a1 * ((Vector3D) u).x;
        y += a1 * ((Vector3D) u).y;
        z += a1 * ((Vector3D) u).z;
    }

    public void PE(IVector u) {
        x += ((Vector3D) u).x;
        y += ((Vector3D) u).y;
        z += ((Vector3D) u).z;
    }

    public void PE(double a) {
        x += a;
        y += a;
        z += a;
    }

    public void ME(IVector u) {
        x -= ((Vector3D) u).x;
        y -= ((Vector3D) u).y;
        z -= ((Vector3D) u).z;
    }

    public void TE(double a) {
        x *= a;
        y *= a;
        z *= a;
    }

    public void TE(IVector u) {
        x *= ((Vector3D) u).x;
        y *= ((Vector3D) u).y;
        z *= ((Vector3D) u).z;
    }

    public void DE(IVector u) {
        x /= ((Vector3D) u).x;
        y /= ((Vector3D) u).y;
        z /= ((Vector3D) u).z;
    }

    public void Ev1Pv2(IVector u1, IVector u2) {
        x = ((Vector3D) u1).x + ((Vector3D) u2).x;
        y = ((Vector3D) u1).y + ((Vector3D) u2).y;
        z = ((Vector3D) u1).z + ((Vector3D) u2).z;
    }

    public void Ev1Mv2(IVector u1, IVector u2) {
        x = ((Vector3D) u1).x - ((Vector3D) u2).x;
        y = ((Vector3D) u1).y - ((Vector3D) u2).y;
        z = ((Vector3D) u1).z - ((Vector3D) u2).z;
    }
    
    public void mod(IVector u) {
        mod((Vector3D) u);
    }

    public void mod(Vector3D u) {
        while (x > u.x)
            x -= u.x;
        while (x < 0.0)
            x += u.x;
        while (y > u.y)
            y -= u.y;
        while (y < 0.0)
            y += u.y;
        while (z > u.z)
            z -= u.z;
        while (z < 0.0)
            z += u.z;
    }

    public double squared() {
        return x * x + y * y + z * z;
    }

    public double Mv1Squared(IVector u) {
        double dx = x - ((Vector3D) u).x;
        double dy = y - ((Vector3D) u).y;
        double dz = z - ((Vector3D) u).z;
        return dx * dx + dy * dy + dz * dz;
    }

    public double dot(IVector u) {
        return x * ((Vector3D) u).x + y * ((Vector3D) u).y + z
                * ((Vector3D) u).z;
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
        x = random.nextDouble() - 0.5;
        y = random.nextDouble() - 0.5;
        z = random.nextDouble() - 0.5;
    }

    public void setRandomInSphere(IRandom random) {//check before using
        double z1 = 0.0;
        double z2 = 0.0;
        double z3 = 0.0;
        double rsq;
        do {

            z1 = 1.0 - 2.0 * Simulation.random.nextDouble();
            z2 = 1.0 - 2.0 * Simulation.random.nextDouble();
            z3 = 1.0 - 2.0 * Simulation.random.nextDouble();

            rsq = z1 * z1 + z2 * z2 + z3 * z3;
        } while (rsq > 1.0);
        x = z1;
        y = z2;
        z = z3;
    }

    /**
     * Creating a random unit vector on unit sphere Uses only two random number
     * generator at a time
     * 
     * Based on M.P. Allen and D.J. Tildesley, Computer Simulation of Liquids, p 349.
     * 
     * @author Jayant Singh
     */
    public void setRandomSphere(IRandom random) {
        double z1, z2, zsq;
        do  {
            z1 = 2.0 * Simulation.random.nextDouble() - 1.0;
            z2 = 2.0 * Simulation.random.nextDouble() - 1.0;
            zsq = z1 * z1 + z2 * z2;
        } while (zsq > 1.0);

        double ranh = 2.0 * Math.sqrt(1.0 - zsq);
        x = z1 * ranh;
        y = z2 * ranh;
        z = 1.0 - 2.0 * zsq;
    }

    public void XE(Vector3D u) {//cross product
        double xNew = y * u.z - z * u.y;
        double yNew = z * u.x - x * u.z;
        z = x * u.y - y * u.x;
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

    public void map(Function function) {
        x = function.f(x);
        y = function.f(y);
        z = function.f(z);
    }
}