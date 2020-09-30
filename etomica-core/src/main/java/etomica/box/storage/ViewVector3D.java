package etomica.box.storage;

import etomica.math.function.IFunction;
import etomica.space.Vector;
import etomica.space3d.Vector3D;
import etomica.util.random.IRandom;

import java.util.Arrays;
import java.util.Objects;
import java.util.StringJoiner;
import java.util.function.DoubleBinaryOperator;

public class ViewVector3D implements Vector {
    private int index;
    private double[] data;

    ViewVector3D(int index, double[] data) {
        this.index = index;
        this.data = data;
    }

    void setIndex(int i) {
        this.index = i;
    }

    void setData(double[] data) {
        this.data = data;
    }

    public double x() {
        return this.data[index];
    }

    public double y() {
        return this.data[index + 1];
    }

    public double z() {
        return this.data[index + 2];
    }

    public void xSet(double d) {
        this.data[index] = d;
    }

    public void ySet(double d) {
        this.data[index + 1] = d;
    }

    public void zSet(double d) {
        this.data[index + 2] = d;
    }

    @Override
    public int getD() {
        return 3;
    }

    @Override
    public void assignTo(double[] array) {
        System.arraycopy(this.data, index, array, 0, 3);
    }

    @Override
    public boolean equals(Vector v) {
        return x() == v.x() && y() == v.y() && z() == v.z();
    }

    @Override
    public double getX(int i) {
        if (i < 0 || i > 2) {
            throw new IllegalArgumentException();
        }
        return this.data[index + i];
    }

    @Override
    public double squared() {
        return this.dot(this);
    }

    @Override
    public boolean isZero() {
        return x() == 0 && y() == 0 && z() == 0;
    }

    @Override
    public double dot(Vector u) {
        return x() * u.x() + y() * u.y() + z() * u.z();
    }

    @Override
    public boolean isNaN() {
        return Double.isNaN(x()) || Double.isNaN(y()) || Double.isNaN(z());
    }

    @Override
    public double Mv1Squared(Vector v1) {
        double dx = x() - v1.x();
        double dy = y() - v1.y();
        double dz = z() - v1.z();
        return dx * dx + dy * dy + dz * dz;
    }

    @Override
    public void setX(int i, double d) {
        if (i < 0 || i > 2) {
            throw new IllegalArgumentException();
        }
        this.data[index + i] = d;
    }

    public void set(double x, double y, double z) {
        xSet(x);
        ySet(y);
        zSet(z);
    }

    @Override
    public void E(Vector u) {
        xSet(u.x());
        ySet(u.y());
        zSet(u.z());
    }

    @Override
    public void E(double a) {
        Arrays.fill(this.data, this.index, this.index + 3, a);
    }

    @Override
    public void E(double... a) {
        System.arraycopy(a, 0, this.data, this.index, 3);
    }

    @Override
    public void PE(Vector u) {
        this.mapWith(u, Double::sum);
    }

    @Override
    public void PE(double a) {
        this.map(x -> x + a);
    }

    @Override
    public void ME(Vector u) {
        this.mapWith(u, (x1, x2) -> x1 - x2);
    }

    @Override
    public void TE(Vector u) {
        this.mapWith(u, (x1, x2) -> x1 * x2);
    }

    @Override
    public void DE(Vector u) {
        this.mapWith(u, (x1, x2) -> x1 / x2);
    }

    @Override
    public void TE(double a) {
        this.map(x -> x * a);
    }

    @Override
    public void Ea1Tv1(double a, Vector v1) {
        this.data[index] = v1.x() * a;
        this.data[index + 1] = v1.y() * a;
        this.data[index + 2] = v1.z() * a;
    }

    @Override
    public void PEa1Tv1(double a, Vector v1) {
        this.mapWith(v1, (x1, x2) -> x1 + (x2 * a));
    }

    @Override
    public void Ev1Pv2(Vector v1, Vector v2) {
        this.data[index] = v1.x() + v2.x();
        this.data[index + 1] = v1.y() + v2.y();
        this.data[index + 2] = v1.z() + v2.z();
    }

    @Override
    public void Ev1Mv2(Vector v1, Vector v2) {
        this.data[index] = v1.x() - v2.x();
        this.data[index + 1] = v1.y() - v2.y();
        this.data[index + 2] = v1.z() - v2.z();
    }

    @Override
    public void mod(Vector u) {
        double x = x();
        double y = y();
        double z = z();
        double ux = u.x();
        double uy = u.y();
        double uz = u.z();
        while (x > ux)
            x -= ux;
        while (x < 0.0)
            x += ux;
        while (y > uy)
            y -= uy;
        while (y < 0.0)
            y += uy;
        while (z > uz)
            z -= uz;
        while (z < 0.0)
            z += uz;
        this.set(x, y, z);
    }

    @Override
    public void normalize() {
        double norm = 1.0 / Math.sqrt(this.squared());
        this.TE(norm);
    }

    @Override
    public void map(IFunction f) {
        this.data[index] = f.f(this.data[index]);
        this.data[index + 1] = f.f(this.data[index + 1]);
        this.data[index + 2] = f.f(this.data[index + 2]);
    }

    @Override
    public void XE(Vector u) {
        set(
                y() * u.z() - z() * u.y(),
                z() * u.x() - x() * u.z(),
                x() * u.y() - y() * u.x()
        );
    }

    /**
     * Creating a random unit vector on unit sphere Uses only two random number
     * generator at a time
     * <p>
     * Based on M.P. Allen and D.J. Tildesley, Computer Simulation of Liquids, p 349.
     * (taken from Marsaglia, G. (1972), Choosing a point from the surface of a sphere,
     * Ann. math. Stat. 43, 645-646.)
     */
    @Override
    public void setRandomSphere(IRandom random) {
        double z1, z2, zsq;
        do {
            z1 = 2.0 * random.nextFixedDouble() - 1.0;
            z2 = 2.0 * random.nextFixedDouble() - 1.0;
            zsq = z1 * z1 + z2 * z2;
        } while (zsq > 1.0);

        double ranh = 2.0 * Math.sqrt(1.0 - zsq);
        data[index] = z1 * ranh;
        data[index + 1] = z2 * ranh;
        data[index + 2] = 1.0 - 2.0 * zsq;
    }

    @Override
    public void setRandomCube(IRandom random) {
        data[index] = random.nextFixedDouble() - 0.5;
        data[index + 1] = random.nextFixedDouble() - 0.5;
        data[index + 2] = random.nextFixedDouble() - 0.5;
    }

    //generate point on surface of sphere according to method (4) here:
    //http://www.math.niu.edu/~rusin/known-math/96/sph.rand  (403)
    //and scale into the interior according to r^2
    @Override
    public void setRandomInSphere(IRandom random) {
        double r = Vector3D.cubeRoot(random.nextFixedDouble());
        double u, v, s;
        do {
            u = 1.0 - 2.0 * random.nextFixedDouble();
            v = 1.0 - 2.0 * random.nextFixedDouble();
            s = u * u + v * v;
        } while (s > 1);
        double ra = 2. * r * Math.sqrt(1. - s);
        data[index] = ra * u;
        data[index + 1] = ra * v;
        data[index + 2] = r * (2 * s - 1.);
    }

    @Override
    public void nearestImage(Vector dimensions) {
        double x = x();
        double y = y();
        double z = z();
        double dx = dimensions.x();
        double dy = dimensions.y();
        double dz = dimensions.z();
        double halfX = dimensions.x() / 2;
        double halfY = dimensions.y() / 2;
        double halfZ = dimensions.z() / 2;
        while (x > halfX)
            x -= dx;
        while (x < -halfX)
            x += dx;
        while (y > halfY)
            y -= dy;
        while (y < -halfY)
            y += dy;
        while (z > halfZ)
            z -= dz;
        while (z < -halfZ)
            z += dz;
        this.set(x, y, z);
    }

    @Override
    public Vector duplicate() {
        return new Vector3D(x(), y(), z());
    }

    public void mapWith(Vector v, DoubleBinaryOperator fn) {
        this.data[index] = fn.applyAsDouble(this.data[index], v.x());
        this.data[index + 1] = fn.applyAsDouble(this.data[index + 1], v.y());
        this.data[index + 2] = fn.applyAsDouble(this.data[index + 2], v.z());
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

    @Override
    public int hashCode() {
        return Objects.hash(x(), y(), z());
    }

    @Override
    public String toString() {
        return new StringJoiner(", ", ViewVector3D.class.getSimpleName() + "(", ")")
                .add(Double.toString(x()))
                .add(Double.toString(y()))
                .add(Double.toString(z()))
                .toString();
    }
}
