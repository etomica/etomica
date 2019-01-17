package etomica.vecarray;

import etomica.math.function.IFunction;
import etomica.space.Vector;
import etomica.util.random.IRandom;

import java.util.Arrays;

public class VecArray {
    private double[] data;
    private VectorFacade[] vectors;
    private int length;
    private static final int D = 3;

    public VecArray(int length) {
        this.data = new double[length * 3];
        this.length = length;

        this.vectors = new VectorFacade[length];
        for (int i = 0; i < length; i++) {
            vectors[i] = new VectorFacade(this, i);
        }
    }

    public int getLength() { return length; }

    public void resize(int newLength) {
        this.length = newLength;
        this.data = Arrays.copyOf(data, newLength * 3);

        // or keep same objects?
        this.vectors = new VectorFacade[newLength];
        for (int i = 0; i < length; i++) {
            vectors[i] = new VectorFacade(this, i * 3);
        }
    }

    public Vector getMut(int i) {
        return vectors[i];
    }

    public void setAll(double x) {
        Arrays.fill(this.data, x);
    }

    public void add(int i, Vector v) {
        data[3 * i] += v.getX(0);
        data[3 * i + 1] += v.getX(1);
        data[3 * i + 2] += v.getX(2);
    }

    public void sub(int i, Vector v) {
        data[3 * i] -= v.getX(0);
        data[3 * i + 1] -= v.getX(1);
        data[3 * i + 2] -= v.getX(2);
    }

    private static class VectorFacade implements Vector {
        private final VecArray array;
        private final double[] data;
        private final int cursor;
        private VectorFacade(VecArray array, int cursor) {
            this.array = array;
            this.cursor = cursor;
            this.data = array.data;
        }

        @Override
        public int getD() {
            return D;
        }

        @Override
        public void assignTo(double[] array) {
            array[0] = this.data[cursor];
            array[1] = this.data[cursor + 1];
            array[2] = this.data[cursor + 2];
        }

        @Override
        public boolean equals(Vector v) {
            return v.getX(0) == data[cursor]
                    && v.getX(1) == data[cursor + 1]
                    && v.getX(2) == data[cursor + 2];
        }

        @Override
        public double getX(int i) {
            return data[cursor + i];
        }

        @Override
        public double squared() {
            return 0;
        }

        @Override
        public boolean isZero() {
            return false;
        }

        @Override
        public double dot(Vector u) {
            return 0;
        }

        @Override
        public boolean isNaN() {
            return false;
        }

        @Override
        public double Mv1Squared(Vector v1) {
            return 0;
        }

        @Override
        public void setX(int i, double d) {

        }

        @Override
        public void E(Vector u) {

        }

        @Override
        public void E(double a) {

        }

        @Override
        public void E(double[] a) {

        }

        @Override
        public void PE(Vector u) {

        }

        @Override
        public void PE(double a) {

        }

        @Override
        public void ME(Vector u) {

        }

        @Override
        public void TE(Vector u) {

        }

        @Override
        public void DE(Vector u) {

        }

        @Override
        public void TE(double a) {

        }

        @Override
        public void Ea1Tv1(double a, Vector v1) {

        }

        @Override
        public void PEa1Tv1(double a, Vector v1) {

        }

        @Override
        public void Ev1Pv2(Vector v1, Vector v2) {

        }

        @Override
        public void Ev1Mv2(Vector v1, Vector v2) {

        }

        @Override
        public void mod(Vector u) {

        }

        @Override
        public void normalize() {

        }

        @Override
        public void map(IFunction f) {

        }

        @Override
        public void XE(Vector u) {

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
            throw new UnsupportedOperationException("This isn't a scratch vector, probably shouldn't be nearest-imaging this");
        }
    }
}
