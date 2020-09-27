package etomica.box.storage;

import etomica.math.function.IFunction;
import etomica.space.Vector;
import etomica.space3d.Vector3D;
import etomica.util.random.IRandom;

import java.util.Arrays;
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
        return x() == v.getX(0) && y() == v.getX(1) && z() == v.getX(2);
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
        return x() * u.getX(0) + y() * u.getX(1) + z() * getX(2);
    }

    @Override
    public boolean isNaN() {
        return Double.isNaN(x()) || Double.isNaN(y()) || Double.isNaN(z());
    }

    @Override
    public double Mv1Squared(Vector v1) {
        double dx = x() - v1.getX(0);
        double dy = x() - v1.getX(1);
        double dz = x() - v1.getX(2);
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
        xSet(u.getX(0));
        ySet(u.getX(1));
        zSet(u.getX(2));
    }

    @Override
    public void E(double a) {
        Arrays.fill(this.data, this.index, this.index + 3, a);
    }

    @Override
    public void E(double[] a) {
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
        this.data[index] = v1.getX(0) * a;
        this.data[index + 1] = v1.getX(1) * a;
        this.data[index + 2] = v1.getX(2) * a;
    }

    @Override
    public void PEa1Tv1(double a, Vector v1) {
        this.mapWith(v1, (x1, x2) -> x1 + (x2 * a));
    }

    @Override
    public void Ev1Pv2(Vector v1, Vector v2) {
        this.data[index] = v1.getX(0) + v2.getX(0);
        this.data[index + 1] = v1.getX(1) + v2.getX(1);
        this.data[index + 2] = v1.getX(2) + v2.getX(2);
    }

    @Override
    public void Ev1Mv2(Vector v1, Vector v2) {
        this.data[index] = v1.getX(0) - v2.getX(0);
        this.data[index + 1] = v1.getX(1) - v2.getX(1);
        this.data[index + 2] = v1.getX(2) - v2.getX(2);
    }

    @Override
    public void mod(Vector u) {
        double x = x();
        double y = y();
        double z = z();
        double ux = u.getX(0);
        double uy = u.getX(1);
        double uz = u.getX(2);
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
                x() * u.getX(1) - y() * u.getX(0),
                y() * u.getX(2) - z() * u.getX(1),
                z() * u.getX(0) - x() * u.getX(2)
        );
    }

    @Override
    public void setRandomSphere(IRandom random) {

    }

    @Override
    public void setRandomCube(IRandom random) {

    }

    @Override
    public void setRandomInSphere(IRandom random) {

    }

    @Override
    public void nearestImage(Vector dimensions) {
        double x = x();
        double y = y();
        double z = z();
        double ux = dimensions.getX(0) / 2;
        double uy = dimensions.getX(1) / 2;
        double uz = dimensions.getX(2) / 2;
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
    public Vector duplicate() {
        return new Vector3D(x(), y(), z());
    }

    public void mapWith(Vector v, DoubleBinaryOperator fn) {
        this.data[index] = fn.applyAsDouble(this.data[index], v.getX(0));
        this.data[index + 1] = fn.applyAsDouble(this.data[index + 1], v.getX(1));
        this.data[index + 2] = fn.applyAsDouble(this.data[index + 2], v.getX(2));
    }
}
