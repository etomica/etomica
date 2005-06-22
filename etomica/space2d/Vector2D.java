package etomica.space2d;

import etomica.Simulation;
import etomica.math.SpecialFunctions;
import etomica.space.Boundary;
import etomica.space.Tensor;
import etomica.space.Vector;
import etomica.space3d.Vector3D;
import etomica.utility.Function;

/*
 * History Created on Jan 24, 2005 by kofke
 */
public final class Vector2D extends etomica.space.Vector { //declared final for

    // efficient method
    // calls

    public static final Vector2D ORIGIN = new Vector2D(0.0, 0.0); //anything
    // using WORK
    // is not
    // thread-safe
    public static final Vector2D WORK = new Vector2D();
    double x, y;

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

    public double[] toArray() {
        return new double[] { x, y };
    }

    public boolean equals(Vector v) {
        return (x == ((Vector2D) v).x) && (y == ((Vector2D) v).y);
    }

    public boolean isZero() {
        return (x == 0.0) && (y == 0.0);
    }

    public void sphericalCoordinates(double[] result) {
        result[0] = Math.sqrt(x * x + y * y);
        result[1] = Math.atan2(y, x); //theta
    }

    public int D() {
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

    public void E(Vector u) {
        x = ((Vector2D) u).x;
        y = ((Vector2D) u).y;
    }

    public void E(double[] u) {
        x = u[0];
        y = u[1];
    } //should check length of array for exception

    public void E(int[] u) {
        x = u[0];
        y = u[1];
    } //should check length of array for exception

    public void E(double a) {
        x = a;
        y = a;
    }

    public void E(double a, double b) {
        x = a;
        y = b;
    }

    public void Ea1Tv1(double a1, Vector u) {
        Vector2D u1 = (Vector2D) u;
        x = a1 * u1.x;
        y = a1 * u1.y;
    }

    public void Ev1Pa1Tv2(Vector v1, double a1, Vector v2) {
        x = ((Vector2D) v1).x + a1 * ((Vector2D) v2).x;
        y = ((Vector2D) v1).y + a1 * ((Vector2D) v2).y;
    }

    public void PEa1Tv1(double a1, Vector u) {
        Vector2D u1 = (Vector2D) u;
        x += a1 * u1.x;
        y += a1 * u1.y;
    }

    public void PE(Vector u) {
        x += ((Vector2D) u).x;
        y += ((Vector2D) u).y;
    }

    public void PE(double a) {
        x += a;
        y += a;
    }

    public void ME(Vector u) {
        x -= ((Vector2D) u).x;
        y -= ((Vector2D) u).y;
    }

    public void PE(int i, double a) {
        if (i == 0)
            x += a;
        else
            y += a;
    }

    public void TE(double a) {
        x *= a;
        y *= a;
    }

    public void TE(Vector u) {
        x *= ((Vector2D) u).x;
        y *= ((Vector2D) u).y;
    }

    public void TE(int i, double a) {
        if (i == 0)
            x *= a;
        else
            y *= a;
    }

    public void DE(double a) {
        x /= a;
        y /= a;
    }

    public void DE(Vector u) {
        x /= ((Vector2D) u).x;
        y /= ((Vector2D) u).y;
    }

    public double Mv1Squared(etomica.space.Vector u) {
        Vector2D u1 = (Vector2D) u;
        double dx = x - u1.x;
        double dy = y - u1.y;
        return dx * dx + dy * dy;
    }

    public void PEa1SGNv1(double a1, Vector v1) {
        x += a1 * SpecialFunctions.sgn(((Vector2D) v1).x);
        y += a1 * SpecialFunctions.sgn(((Vector2D) v1).y);
    }

    public void Ev1Pv2(Vector u1, Vector u2) {
        Vector2D v1 = (Vector2D) u1;
        Vector2D v2 = (Vector2D) u2;
        x = v1.x + v2.x;
        y = v1.y + v2.y;
    }

    public void Ev1Mv2(Vector u1, Vector u2) {
        Vector2D v1 = (Vector2D) u1;
        Vector2D v2 = (Vector2D) u2;
        x = v1.x - v2.x;
        y = v1.y - v2.y;
    }

    public void mod(Vector u) {
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

    public void mod(double a) {
        while (x > a)
            x -= a;
        while (x < 0.0)
            x += a;
        while (y > a)
            y -= a;
        while (y < 0.0)
            y += a;
    }

    //		public void EModShift(Space.Vector r, Space.Vector u) {
    //			EModShift((Vector)r, (Vector)u);
    //		}
    //sets this equal to (r mod u) - r
    public void EModShift(Vector r, Vector u) {
        Vector2D r2d = (Vector2D) r;
        Vector2D u2d = (Vector2D) u;
        x = r2d.x;
        while (x > u2d.x)
            x -= u2d.x;
        while (x < 0.0)
            x += u2d.x;
        x -= r2d.x;
        y = r2d.y;
        while (y > u2d.y)
            y -= u2d.y;
        while (y < 0.0)
            y += u2d.y;
        y -= r2d.y;
    }

    public void EMod2Shift(Vector r, Vector u) {
        Vector2D r2d = (Vector2D) r;
        Vector2D u2d = (Vector2D) u;
        x = r2d.x;
        while (x > +u2d.x)
            x -= (u2d.x + u2d.x);
        while (x < -u2d.x)
            x += (u2d.x + u2d.x);
        x -= r2d.x;
        y = r2d.y;
        while (y > +u2d.y)
            y -= (u2d.y + u2d.y);
        while (y < -u2d.y)
            y += (u2d.y + u2d.y);
        y -= r2d.y;
    }

    public etomica.space.Vector P(Vector u) {
        Vector2D u1 = (Vector2D) u;
        WORK.x = x + u1.x;
        WORK.y = y + u1.y;
        return WORK;
    }

    public etomica.space.Vector M(Vector u) {
        Vector2D u1 = (Vector2D) u;
        WORK.x = x - u1.x;
        WORK.y = y - u1.y;
        return WORK;
    }

    public etomica.space.Vector T(Vector u) {
        Vector2D u1 = (Vector2D) u;
        WORK.x = x * u1.x;
        WORK.y = y * u1.y;
        return WORK;
    }

    public etomica.space.Vector D(Vector u) {
        Vector2D u1 = (Vector2D) u;
        WORK.x = x / u1.x;
        WORK.y = y / u1.y;
        return WORK;
    }

    public void abs() {
        x = (x > 0) ? x : -x;
        y = (y > 0) ? y : -y;
    }

    public double min() {
        return (x < y) ? x : y;
    }

    public double max() {
        return (x > y) ? x : y;
    }

    public double squared() {
        return x * x + y * y;
    }

    public double dot(Vector u) {
        return x * ((Vector2D) u).x + y * ((Vector2D) u).y;
    }

    public Vector3D cross(Vector3D u) {
        return new Vector3D(y * u.x(2), -x * u.x(2), x * u.x(1) - y * u.x(0));
        //        	work.setX(0, y*u.x(2));
        //            work.setX(1,-x*u.x(2));
        //            work.setX(2, x*u.x(1) - y*u.x(0));
    }

    public Vector3D cross(Vector2D u) {
        return new Vector3D(0.0, 0.0, x * u.y - y * u.x);
        //            Space3D.Vector.WORK.setX(0, 0.0);
        //            Space3D.Vector.WORK.setX(1, 0.0);
        //            Space3D.Vector.WORK.setX(2, x*u.y - y*u.x);
        //            return Space3D.Vector.WORK;
    }

    /**
     * Replaces this vector with its cross-product with the given 3D vector,
     * with result projected onto the 2D plane. This vector becomes the result
     * of (this vector) X u.
     */
    public void XE(Vector3D u) {
        double xNew = y * u.x(2);
        y = -x * u.x(2);
        x = xNew;
    }

    public void normalize() {
        double norm = Math.sqrt(1 / (x * x + y * y));
        x *= norm;
        y *= norm;
    }

    public void transform(Tensor A) {
        double x0 = x;
        double y0 = y;
        x = ((Tensor2D) A).xx * x0 + ((Tensor2D) A).xy * y0;
        y = ((Tensor2D) A).yx * x0 + ((Tensor2D) A).yy * y0;
    }

    public void transform(Boundary b, Vector r, Tensor A) {
        Vector2D r2 = (Vector2D) r;
        Tensor2D A2 = (Tensor2D) A;
        WORK.x = x - r2.x;
        WORK.y = y - r2.y;
        b.nearestImage(WORK);
        x = r2.x + A2.xx * WORK.x + A2.xy * WORK.y;
        y = r2.y + A2.yx * WORK.x + A2.yy * WORK.y;
    }

    public void randomStep(double d) {
        x += (2. * Simulation.random.nextDouble() - 1.0) * d;
        y += (2. * Simulation.random.nextDouble() - 1.0) * d;
    } //uniformly distributed random step in x and y, within +/- d

    public void setRandom(double d) {
        x = Simulation.random.nextDouble() * d;
        y = Simulation.random.nextDouble() * d;
    }

    public void setRandom(double dx, double dy) {
        x = Simulation.random.nextDouble() * dx;
        y = Simulation.random.nextDouble() * dy;
    }

    public void setRandom(Vector u) {
        setRandom(((Vector2D) u).x, ((Vector2D) u).y);
    }

    public void setRandomCube() {
        x = Simulation.random.nextDouble() - 0.5;
        y = Simulation.random.nextDouble() - 0.5;
    }

    public void setRandomSphere() {
        x = Math.cos(2 * Math.PI * Simulation.random.nextDouble());
        y = Math.sqrt(1.0 - x * x);
        if (Simulation.random.nextDouble() < 0.5)
            y = -y;
    }

    // random point in a unit sphere
    public void setRandomInSphere() {//check before using
        double z1 = 0.0;
        double z2 = 0.0;
        double rsq;
        do {
            z1 = 1.0 - 2.0 * Simulation.random.nextDouble();
            z2 = 1.0 - 2.0 * Simulation.random.nextDouble();
            rsq = z1 * z1 + z2 * z2;
        } while (rsq > 1.0);
        x = z1;
        y = z2;
    }

    public void randomRotate(double thetaStep) {
        double deltheta = (2 * Simulation.random.nextDouble() - 1.0)
                * thetaStep;
        double theta = Math.atan2(y, x);
        theta += deltheta;
        double r = Math.sqrt(x * x + y * y);
        x = r * Math.cos(theta);
        y = r * Math.sin(theta);
    }

    public boolean isNaN() {
        return Double.isNaN(x) || Double.isNaN(y);
    }

    public void map(Function function) {
        x = function.f(x);
        y = function.f(y);
    }

}