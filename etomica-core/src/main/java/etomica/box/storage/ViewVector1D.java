package etomica.box.storage;

import etomica.math.function.IFunction;
import etomica.space.Vector;
import etomica.space2d.Vector2D;
import etomica.util.random.IRandom;

import java.util.Objects;
import java.util.StringJoiner;
import java.util.function.DoubleBinaryOperator;

public class ViewVector1D implements Vector {
    private int index;
    private double[] data;

    ViewVector1D(int index, double[] data) {
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
        throw new RuntimeException("Not a 2D vector");
    }

    public double z() {
        throw new RuntimeException("Not a 3D vector");
    }

    public void xSet(double d) {
        this.data[index] = d;
    }

    public void ySet(double d) {
        throw new RuntimeException("Not a 2D vector");
    }

    public void zSet(double d) {
        throw new RuntimeException("Not a 3D vector");
    }

    @Override
    public int getD() {
        return 1;
    }

    @Override
    public void assignTo(double[] array) {
        array[0] = data[index];
    }

    @Override
    public boolean equals(Vector v) {
        return x() == v.x();
    }

    @Override
    public double getX(int i) {
        if (i != 0) {
            throw new IllegalArgumentException();
        }
        return this.data[index];
    }

    @Override
    public double squared() {
        return this.dot(this);
    }

    @Override
    public boolean isZero() {
        return x() == 0;
    }

    @Override
    public double dot(Vector u) {
        return x() * u.x();
    }

    @Override
    public boolean isNaN() {
        return Double.isNaN(x());
    }

    @Override
    public double Mv1Squared(Vector v1) {
        double dx = x() - v1.x();
        return dx * dx;
    }

    @Override
    public void setX(int i, double d) {
        if (i != 0) {
            throw new IllegalArgumentException();
        }
        this.data[index] = d;
    }

    @Override
    public void E(Vector u) {
        xSet(u.x());
    }

    @Override
    public void E(double a) {
        data[index] = a;
    }

    @Override
    public void E(double... a) {
        data[index] = a[0];
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
    }

    @Override
    public void PEa1Tv1(double a, Vector v1) {
        this.mapWith(v1, (x1, x2) -> x1 + (x2 * a));
    }

    @Override
    public void Ev1Pv2(Vector v1, Vector v2) {
        this.data[index] = v1.x() + v2.x();
    }

    @Override
    public void Ev1Mv2(Vector v1, Vector v2) {
        this.data[index] = v1.x() - v2.x();
    }

    @Override
    public void mod(Vector u) {
        double x = x();
        double ux = u.x();
        while (x > ux)
            x -= ux;
        while (x < 0.0)
            x += ux;
        data[index] = x;
    }

    @Override
    public void normalize() {
        double norm = 1.0 / Math.sqrt(this.squared());
        this.TE(norm);
    }

    @Override
    public void map(IFunction f) {
        this.data[index] = f.f(this.data[index]);
    }

    @Override
    public void XE(Vector u) {
        throw new RuntimeException("Cross product only exists for 3D vectors");
    }

    @Override
    public void setRandomSphere(IRandom random) {
        data[index] = random.nextInt(2) == 0 ? -1 : +1;
    }

    @Override
    public void setRandomCube(IRandom random) {
        data[index] = random.nextFixedDouble() - 0.5;
    }

    @Override
    public void setRandomInSphere(IRandom random) {
        data[index] = 1.0 - 2.0 * random.nextFixedDouble();
    }

    @Override
    public void nearestImage(Vector dimensions) {
        double x = x();
        double dx = dimensions.x();
        double halfX = dimensions.x() / 2;
        while (x > halfX)
            x -= dx;
        while (x < -halfX)
            x += dx;
        data[index] = x;
    }

    @Override
    public Vector duplicate() {
        return new Vector2D(x(), y());
    }

    public void mapWith(Vector v, DoubleBinaryOperator fn) {
        this.data[index] = fn.applyAsDouble(this.data[index], v.x());
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null) return false;
        if (!(o instanceof Vector) || ((Vector) o).getD() != this.getD()) {
            return false;
        }
        Vector v = ((Vector) o);
        return v.x() == this.x();
    }

    @Override
    public int hashCode() {
        return Objects.hash(x());
    }

    @Override
    public String toString() {
        return new StringJoiner(", ", ViewVector1D.class.getSimpleName() + "(", ")")
                .add(Double.toString(x()))
                .toString();
    }
}
