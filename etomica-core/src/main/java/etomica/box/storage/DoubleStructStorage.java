package etomica.box.storage;

import java.lang.reflect.Array;
import java.util.Arrays;

public abstract class DoubleStructStorage<V> implements Storage {
    protected double[] data;
    private V[] views;
    private int count;
    protected final int stride;

    private static final double SIZE_INCREASE_RATIO = 0.3;

    public DoubleStructStorage(int stride, int count, Class<? extends V> viewClass) {
        this.stride = stride;
        this.count = count;
        this.data = new double[count * stride];

        @SuppressWarnings("unchecked")
        final V[] views = (V[]) Array.newInstance(viewClass, count);
        this.views = views;
    }

    protected abstract V makeView(int i);

    protected abstract void updateIndex(V view, int newIdx);

    protected abstract void updateData(V[] views, double[] newData);

    public V get(int i) {
        return views[i];
    }

    public V create(int i) {
        if (this.views[i] != null) {
            throw new IllegalArgumentException("Already exists: " + i);
        }
        V view = makeView(i);
        this.views[i] = view;
        return view;
    }

    private void resize(int newSize) {
        this.data = Arrays.copyOf(this.data, newSize * stride);
        this.views = Arrays.copyOf(this.views, newSize);
        this.updateData(views, data);
    }

    public void addNull(int n) {
        if (this.count + n > this.views.length) {
            resize(Math.max((int) (count * (1.0 + SIZE_INCREASE_RATIO) + 1), count + n));
        }
        count += n;
    }

    @Override
    public void add(int n) {
        addNull(n);
        for (int i = count - n; i < count; i++) {
            this.views[i] = makeView(i);
        }
    }

    @Override
    public void ensureCapacity(int expectedAdditions) {
        if (count + expectedAdditions > views.length) {
            resize(count + expectedAdditions);
        }
    }

    @Override
    public void swapRemove(int i) {
        if (views[count - 1] != null) {
            if (i != count - 1) {
                System.arraycopy(data, (count - 1) * stride, data, i * stride, stride);
                views[i] = views[count - 1];
                this.updateIndex(views[i], i);
            }
            this.views[count - 1] = null;
        } else {
            views[i] = null;
        }
        this.count--;
    }

    @Override
    public int size() {
        return this.count;
    }
}
