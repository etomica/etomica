package etomica.box.storage;

import etomica.math.function.IFunction;
import etomica.space.Vector;
import etomica.space2d.Vector2D;
import etomica.util.random.IRandom;

import java.util.Arrays;
import java.util.Objects;
import java.util.StringJoiner;
import java.util.function.DoubleBinaryOperator;

public class ViewVector2D implements Vector {
    private int index;
    private double[] data;

    ViewVector2D(int index, double[] data) {
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
        throw new RuntimeException("Not a 3D vector");
    }

    public void xSet(double d) {
        this.data[index] = d;
    }

    public void ySet(double d) {
        this.data[index + 1] = d;
    }

    public void zSet(double d) {
        throw new RuntimeException("Not a 3D vector");
    }

    @Override
    public int getD() {
        return 2;
    }

    @Override
    public void assignTo(double[] array) {
        System.arraycopy(this.data, index, array, 0, 2);
    }

    @Override
    public boolean equals(Vector v) {
        return x() == v.x() && y() == v.y();
    }

    @Override
    public double getX(int i) {
        if (i < 0 || i > 1) {
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
        return x() == 0 && y() == 0;
    }

    @Override
    public double dot(Vector u) {
        return x() * u.x() + y() * u.y();
    }

    @Override
    public boolean isNaN() {
        return Double.isNaN(x()) || Double.isNaN(y());
    }

    @Override
    public double Mv1Squared(Vector v1) {
        double dx = x() - v1.x();
        double dy = y() - v1.y();
        return dx * dx + dy * dy;
    }

    @Override
    public void setX(int i, double d) {
        if (i < 0 || i > 1) {
            throw new IllegalArgumentException();
        }
        this.data[index + i] = d;
    }

    public void set(double x, double y) {
        xSet(x);
        ySet(y);
    }

    @Override
    public void E(Vector u) {
        xSet(u.x());
        ySet(u.y());
    }

    @Override
    public void E(double a) {
        Arrays.fill(this.data, this.index, this.index + 2, a);
    }

    @Override
    public void E(double... a) {
        System.arraycopy(a, 0, this.data, this.index, 2);
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
    }

    @Override
    public void PEa1Tv1(double a, Vector v1) {
        this.mapWith(v1, (x1, x2) -> x1 + (x2 * a));
    }

    @Override
    public void Ev1Pv2(Vector v1, Vector v2) {
        this.data[index] = v1.x() + v2.x();
        this.data[index + 1] = v1.y() + v2.y();
    }

    @Override
    public void Ev1Mv2(Vector v1, Vector v2) {
        this.data[index] = v1.x() - v2.x();
        this.data[index + 1] = v1.y() - v2.y();
    }

    @Override
    public void mod(Vector u) {
        double x = x();
        double y = y();
        double ux = u.x();
        double uy = u.y();
        while (x > ux)
            x -= ux;
        while (x < 0.0)
            x += ux;
        while (y > uy)
            y -= uy;
        while (y < 0.0)
            y += uy;
        this.set(x, y);
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
    }

    @Override
    public void XE(Vector u) {
        throw new RuntimeException("Cross product only exists for 3D vectors");
    }

    @Override
    public void setRandomSphere(IRandom random) {
        double theta = 2 * Math.PI * random.nextFixedDouble();
        double x = Math.cos(theta);
        double y = Math.sqrt(1.0 - x * x);
        if (theta > Math.PI) y = -y;
        data[index] = x;
        data[index + 1] = y;
    }

    @Override
    public void setRandomCube(IRandom random) {
        data[index] = random.nextFixedDouble() - 0.5;
        data[index + 1] = random.nextFixedDouble() - 0.5;
    }

    @Override
    public void setRandomInSphere(IRandom random) {
        double z1 = 0.0;
        double z2 = 0.0;
        double rsq;
        do {
            z1 = 1.0 - 2.0 * random.nextFixedDouble();
            z2 = 1.0 - 2.0 * random.nextFixedDouble();
            rsq = z1 * z1 + z2 * z2;
        } while (rsq > 1.0);
        data[index] = z1;
        data[index + 1] = z2;
    }

    @Override
    public void nearestImage(Vector dimensions) {
        double x = x();
        double y = y();
        double dx = dimensions.x();
        double dy = dimensions.y();
        double halfX = dimensions.x() / 2;
        double halfY = dimensions.y() / 2;
        while (x > halfX)
            x -= dx;
        while (x < -halfX)
            x += dx;
        while (y > halfY)
            y -= dy;
        while (y < -halfY)
            y += dy;
        this.set(x, y);
    }

    @Override
    public Vector duplicate() {
        return new Vector2D(x(), y());
    }

    public void mapWith(Vector v, DoubleBinaryOperator fn) {
        this.data[index] = fn.applyAsDouble(this.data[index], v.x());
        this.data[index + 1] = fn.applyAsDouble(this.data[index + 1], v.y());
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null) return false;
        if (!(o instanceof Vector) || ((Vector) o).getD() != this.getD()) {
            return false;
        }
        Vector v = ((Vector) o);
        return v.x() == this.x() && v.y() == this.y();
    }

    @Override
    public int hashCode() {
        return Objects.hash(x(), y());
    }

    @Override
    public String toString() {
        return new StringJoiner(", ", ViewVector2D.class.getSimpleName() + "(", ")")
                .add(Double.toString(x()))
                .add(Double.toString(y()))
                .toString();
    }
}
