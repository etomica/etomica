package etomica.box.storage;

import java.util.Arrays;
import java.util.BitSet;

public class DoubleStorage implements Storage {
    private static final double SIZE_INCREASE_RATIO = 0.3;
    private double[] data;
    private DoubleWrapper[] views;
    private final BitSet validBits;
    private int count;

    public DoubleStorage(int count) {
        this.data = new double[count];
        this.views = new DoubleWrapper[count];
        this.validBits = new BitSet(count);
        this.count = count;
    }

    public double get(int i) {
        return this.data[i];
    }

    public DoubleWrapper create(int i) {
        DoubleWrapper w = new DoubleWrapper();
        w.data = data;
        w.idx = i;
        views[i] = w;
        return w;
    }

    public void set(int i, double d) {
        this.data[i] = d;
    }

    private void resize(int newSize) {
        this.data = Arrays.copyOf(data, newSize);
        this.views = Arrays.copyOf(views, newSize);
        for (DoubleWrapper view : views) {
            if (view != null) {
                view.data = data;
            }
        }
    }

    public int size() {
        return this.data.length;
    }

    public void add(int n) {
        this.addNull(n);
        validBits.set(count - n, count);
    }

    public void addNull(int n) {
        if (this.count + n >= this.data.length) {
            resize(Math.max((int) (count * (1.0 + SIZE_INCREASE_RATIO) + 1), count + n));
        }
        validBits.clear(count, count + n);
        count += n;
    }

    @Override
    public void ensureCapacity(int expectedAdditions) {
        if (count + expectedAdditions > data.length) {
            resize(count + expectedAdditions);
        }
    }

    public void swapRemove(int i) {
        this.data[i] = this.data[count - 1];
        this.views[i] = this.views[count - 1];
        this.views[i].idx = i;
        this.validBits.set(i, validBits.get(count - 1));
        this.validBits.clear(count - 1);
        count--;
    }

    public static class DoubleWrapper {
        private double[] data;
        private int idx;

        public double get() {
            return data[idx];
        }

        public void set(double d) {
            data[idx] = d;
        }
    }
}
