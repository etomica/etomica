package etomica.box.storage;

import java.util.Arrays;

public class DoubleStorage {
    private double[] data;

    public DoubleStorage(int count) {
        this.data = new double[count];
    }

    public double get(int i) {
        return this.data[i];
    }

    public void resize(int newSize) {
        this.data = Arrays.copyOf(data, newSize);
    }

    public int size() {
        return this.data.length;
    }
}
