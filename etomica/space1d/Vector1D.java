package etomica.space1d;

import etomica.math.SpecialFunctions;
import etomica.simulation.Simulation;
import etomica.space.IVector;
import etomica.util.Function;

/**
 * Implementation of the Vector class for a 1-dimensional space. In this case the vector
 * is a trivial object formed from just one element.
 */
public final class Vector1D implements IVector, java.io.Serializable { 
    
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

    public int D() {
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

    public double[] toArray() {
        return new double[] { x };
    }

    public void sphericalCoordinates(double[] result) {
        result[0] = x;
    }

    public void E(double a) {
        x = a;
    }

    public void E(double[] a) {
        x = a[0];
    }

    public void E(int[] a) {
        x = a[0];
    }

    public void Ea1Tv1(double a1, IVector u) {
        Vector1D u1 = (Vector1D) u;
        x = a1 * u1.x;
    }

    public void Ev1Pa1Tv2(IVector v1, double a1, IVector v2) {
        x = ((Vector1D) v1).x + a1 * ((Vector1D) v2).x;
    }

    public void PEa1Tv1(double a1, IVector u) {
        x += a1 * ((Vector1D) u).x;
    }

    public void PE(double a) {
        x += a;
    }

    public void PE(int i, double a) {
        x += a;
    }

    public void TE(double a) {
        x *= a;
    }

    public void TE(int i, double a) {
        x *= a;
    }

    public void DE(double a) {
        x /= a;
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
    
    public void truncate(double eps) {
        if(x < eps && -x < eps) x = 0.0;
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

    public void mod(double a) {
        while (x > a)
            x -= a;
        while (x < 0.0)
            x += a;
    }

    //        public void EModShift(Space.Vector r, Space.Vector u) {
    //        	EModShift((Vector)r, (Vector)u);
    //        }
    //sets this equal to (r mod u) - r
    public void EModShift(IVector r, IVector u) {
        double rx = ((Vector1D) r).x;
        double ux = ((Vector1D) u).x;
        x = rx;
        while (x >= ux)
            x -= ux;
        while (x < 0.)
            x += ux;
        x -= rx;
    }

    public void EMod2Shift(IVector r, IVector u) {
        double rx = ((Vector1D) r).x;
        double ux = ((Vector1D) u).x;
        x = rx;
        while (x > +ux)
            x -= (ux + ux);
        while (x < -ux)
            x += (ux + ux);
        x -= rx;
    }

    public void PEa1SGNv1(double a1, IVector v1) {
        x += a1 * SpecialFunctions.sgn(((Vector1D) v1).x);
    }

    public IVector P(IVector u) {
        Vector1D work = new Vector1D();
        work.x = x + ((Vector1D)u).x;
        return work;
    }

    public IVector M(IVector u) {
        Vector1D work = new Vector1D();
        work.x = x - ((Vector1D)u).x;
        return work;
    }

    public IVector T(IVector u) {
        Vector1D work = new Vector1D();
        work.x = x * ((Vector1D)u).x;
        return work;
    }

    public IVector D(IVector u) {
        Vector1D work = new Vector1D();
        work.x = x / ((Vector1D)u).x;
        return work;
    }

    public void abs() {
        x = (x > 0) ? x : -x;
    }

    public void minE(IVector v) {
        if(((Vector1D)v).x < x) x = ((Vector1D)v).x;
    }

    public void maxE(IVector v) {
        if(((Vector1D)v).x > x) x = ((Vector1D)v).x;
    }

   public double min() {
        return x;
    }

    public double max() {
        return x;
    }

    public double squared() {
        return x * x;
    }

    public void normalize() {
        x = 1.0;
    }

    public void transform(etomica.space.Tensor A) {
        x = ((Tensor1D) A).xx * x;
    }

    public void randomStep(double d) {
        x += (2. * Simulation.random.nextDouble() - 1.0) * d;
    } //uniformly distributed random step in x and y, within +/- d

    public void setRandom(double d) {
        x = Simulation.random.nextDouble() * d;
    }

    public void setRandom(IVector u) {
        setRandom(((Vector1D) u).x);
    }

    public void setRandomCube() {
        x = Simulation.random.nextDouble() - 0.5;
    }

    public void setRandomInSphere() {
        setRandomCube();
    }

    public void setRandomSphere() {
        randomDirection();
    }

    public void randomDirection() {
        x = (Simulation.random.nextBoolean()) ? -1.0 : +1.0;
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

    public void randomRotate(double deltheta) {//no implementation in 1D
    }

    public boolean isNaN() {
        return Double.isNaN(x);
    }

    public void map(Function function) {
        x = function.f(x);
    }

}