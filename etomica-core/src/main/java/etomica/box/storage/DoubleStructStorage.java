package etomica.box.storage;

import java.lang.reflect.Array;
import java.util.Arrays;

public abstract class DoubleStructStorage<V> extends AbstractObjectStorage<V> {
    protected double[] data;
    protected final int stride;

    public DoubleStructStorage(int stride, Class<? extends V> viewClass) {
        super(viewClass);
        this.stride = stride;
        this.data = new double[count * stride];
    }

    protected abstract void updateIndex(V view, int newIdx);

    protected abstract void updateData(V[] views, double[] newData);

    protected void resize(int newSize) {
        super.resize(newSize);
        this.data = Arrays.copyOf(this.data, newSize * stride);
        this.updateData(objects, data);
    }

    @Override
    public void swapRemove(int i) {
        if (objects[count - 1] != null) {
            if (i != count - 1) {
                System.arraycopy(data, (count - 1) * stride, data, i * stride, stride);
                objects[i] = objects[count - 1];
                this.updateIndex(objects[i], i);
            }
            this.objects[count - 1] = null;
        } else {
            objects[i] = null;
        }
        this.count--;
    }
}
