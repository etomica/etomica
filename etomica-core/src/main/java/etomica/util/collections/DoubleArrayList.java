package etomica.util.collections;

import etomica.util.Debug;

import java.util.Arrays;

public final class DoubleArrayList {
    private double[] data;
    private int size;

    public DoubleArrayList(int initialCapacity) {
        super();
        this.data = new double[initialCapacity];
        this.size = 0;
    }

    private void grow() {
        this.data = Arrays.copyOf(this.data, this.data.length * 2);
    }

    public DoubleArrayList() {
        this(32);
    }

    public void ensureCapacity(int capacity) {
        if (this.data.length < capacity) {
            this.data = Arrays.copyOf(this.data, capacity);
        }
    }

    public double getDouble(int i) {
        if (Debug.ON && i >= size) {
            throw new IndexOutOfBoundsException();
        }

        return data[i];
    }

    public void add(double x) {
        if (size == data.length) {
            this.grow();
        }

        data[size] = x;
        size++;
    }

    public void plusEquals(int i, double x) {
        data[i] += x;
    }

    public void replace(int i, double x) {
        if (Debug.ON && i >= size) {
            throw new IndexOutOfBoundsException();
        }
        data[i] = x;
    }

    public void clear() {
        this.size = 0;
    }

    public void setAll(double x) {
        Arrays.fill(data, x);
    }

    public void set(int i, int x) {
        if (Debug.ON && i >= size) throw new IndexOutOfBoundsException();
        data[i] = x;
    }

    public int size() {
        return size;
    }
}

