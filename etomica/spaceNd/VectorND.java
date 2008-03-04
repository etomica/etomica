package etomica.spaceNd;

import etomica.api.IFunction;
import etomica.api.IRandom;
import etomica.api.IVector;
import etomica.space.IVectorRandom;

/**
 * Implementation of the Vector class for a 3-dimensional space.
 */
public final class VectorND implements IVectorRandom, java.io.Serializable {

    protected final double[] x;
    private static final long serialVersionUID = 1L;
    
    public VectorND(int D) {
        x = new double[D];
    }
    
    public int getD() {
        return x.length;
    }

    public VectorND(double[] a) {
        x = new double[a.length];
        for (int i=0; i<a.length; i++) {
            x[i] = a[i];
        }
    }

    public String toString() {
        String s = "(";
        for (int i=0; i<x.length-1; i++) {
            s += x[i] + ",";
        }
        s += x[x.length-1]+ ")";
        return s;
    }

    public double x(int i) {
        return x[i];
    }

    public void assignTo(double[] array) {
        for (int i=0; i<x.length; i++) {
            array[i] = x[i];
        }
    }

    public boolean equals(IVector v) {
        for (int i=0; i<x.length; i++) {
            if (v.x(i) != x[i]) return false;
        }
        return true;
    }

    public boolean isZero() {
        for (int i=0; i<x.length; i++) {
            if (x[i] != 0) return false;
        }
        return true;
    }

    public void E(IVector u) {
        for (int i=0; i<x.length; i++) {
            x[i] = u.x(i);
        }
    }

    public void E(double a) {
        for (int i=0; i<x.length; i++) {
            x[i] = a;
        }
    }

    public void E(double[] u) {
        for (int i=0; i<x.length; i++) {
            x[i] = u[i];
        }
    }
    
    public void Ea1Tv1(double a1, IVector u) {
        for (int i=0; i<x.length; i++) {
            x[i] = a1 * u.x(i);
        }
    }

    public void PEa1Tv1(double a1, IVector u) {
        for (int i=0; i<x.length; i++) {
            x[i] += a1 * u.x(i);
        }
    }

    public void PE(IVector u) {
        for (int i=0; i<x.length; i++) {
            x[i] += u.x(i);
        }
    }

    public void PE(double a) {
        for (int i=0; i<x.length; i++) {
            x[i] += a;
        }
    }

    public void ME(IVector u) {
        for (int i=0; i<x.length; i++) {
            x[i] -= u.x(i);
        }
    }

    public void TE(double a) {
        for (int i=0; i<x.length; i++) {
            x[i] *= a;
        }
    }

    public void TE(IVector u) {
        for (int i=0; i<x.length; i++) {
            x[i] *= u.x(i);
        }
    }

    public void DE(IVector u) {
        for (int i=0; i<x.length; i++) {
            x[i] /= u.x(i);
        }
    }

    public void Ev1Pv2(IVector u1, IVector u2) {
        for (int i=0; i<x.length; i++) {
            x[i] = u1.x(i) + u2.x(i);
        }
    }

    public void Ev1Mv2(IVector u1, IVector u2) {
        for (int i=0; i<x.length; i++) {
            x[i] = u1.x(i) - u2.x(i);
        }
    }
    
    public void mod(final IVector u) {
        for (int i=0; i<x.length; i++) {
            while (x[i] > u.x(i)) {
                x[i] -= u.x(i);
            }
            while (x[i] < 0) {
                x[i] += u.x(i);
            }
        }
    }

    public double squared() {
        double s = 0;
        for (int i=0; i<x.length; i++) {
            s += x[i] * x[i];
        }
        return s;
    }

    public double Mv1Squared(IVector u) {
        double s = 0;
        for (int i=0; i<x.length; i++) {
            double dx = x[i] - u.x(i);
            s += dx * dx;
        }
        return s;
    }

    public double dot(IVector u) {
        double s = 0;
        for (int i=0; i<x.length; i++) {
            s += x[i] * u.x(i);
        }
        return s;
    }

    public void setX(int a, double d) {
        x[a] = d;
    }

    public void setRandomCube(IRandom random) {
        for (int i=0; i<x.length; i++) {
            x[i] = random.nextDouble() - 0.5;
        }
    }

    public void setRandomInSphere(IRandom random) {//check before using
        double rsq;
        do {
            rsq = 0;
            for (int i=0; i<x.length; i++) {
                x[i] = 1.0 - 2.0 * random.nextDouble();
                rsq += x[i] * x[i];
            }
        } while (rsq > 1.0);
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
        if (x.length == 3) {
            double z1, z2, zsq;
            do  {
                z1 = 2.0 * random.nextDouble() - 1.0;
                z2 = 2.0 * random.nextDouble() - 1.0;
                zsq = z1 * z1 + z2 * z2;
            } while (zsq > 1.0);
    
            double ranh = 2.0 * Math.sqrt(1.0 - zsq);
            x[0] = z1 * ranh;
            x[1] = z2 * ranh;
            x[2] = 1.0 - 2.0 * zsq;
        }
        else if (x.length == 2) {
            double theta = random.nextDouble() * Math.PI * 2;
            x[0] = Math.cos(theta);
            x[1] = Math.sin(theta);
        }
        else if (x.length == 1) {
            x[0] = random.nextInt(2) * 2 - 1;
        }
        else {
            throw new RuntimeException("Can only do random sphere for 1, 2, 3-D vectors");
        }
    }

    public void XE(IVector u) {//cross product
        if (x.length == 3) {
            double xNew = x[1] * u.x(2) - x[2] * u.x(1);
            double yNew = x[2] * u.x(0) - x[0] * u.x(2);
            x[2] = x[0] * u.x(2) - x[1] * u.x(0);
            x[1] = yNew;
            x[0] = xNew;
        }
        else {
            throw new RuntimeException("Can only do cross product for 3-D vectors");
        }
    }

    public void normalize() {
        TE(1. / Math.sqrt(squared()));
    }

    public boolean isNaN() {
        for (int i=0; i<x.length; i++) {
            if (Double.isNaN(x[i])) {
                return true;
            }
        }
        return false;
    }

    public void map(IFunction function) {
        for (int i=0; i<x.length; i++) {
            x[i] = function.f(x[i]);
        }
    }
}