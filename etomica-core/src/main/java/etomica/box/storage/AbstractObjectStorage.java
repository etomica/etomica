package etomica.box.storage;

import java.lang.reflect.Array;
import java.util.Arrays;
import java.util.Objects;

public abstract class AbstractObjectStorage<V> implements Storage {
    protected V[] objects;
    protected int count;

    private static final double SIZE_INCREASE_RATIO = 0.3;

    public AbstractObjectStorage(Class<? extends V> cls) {
        this.count = 0;

        @SuppressWarnings("unchecked")
        final V[] data = (V[]) Array.newInstance(cls, count);
        this.objects = data;
    }

    protected abstract V createObject(int idx);

    protected void destroy(V value, int idx) {}

    public V get(int i) {
        return objects[i];
    }

    public V create(int i) {
        if (this.objects[i] != null) {
            throw new IllegalArgumentException("Already exists: " + i);
        }
        if (i >= count) {
            throw new IndexOutOfBoundsException();
        }
        V obj = this.createObject(i);
        this.objects[i] = obj;
        return obj;
    }

    protected void resize(int newSize) {
        this.objects = Arrays.copyOf(this.objects, newSize);
    }

    @Override
    public void add(int n) {
        addNull(n);
        for (int i = count - n; i < count; i++) {
            this.create(i);
        }
    }

    @Override
    public void addNull(int n) {
        if (this.count + n > this.objects.length) {
            resize(Math.max((int) (count * (1.0 + SIZE_INCREASE_RATIO) + 1), count + n));
        }
        count += n;
    }

    @Override
    public void ensureCapacity(int expectedAdditions) {
        if (count + expectedAdditions > objects.length) {
            resize(count + expectedAdditions);
        }
    }

    @Override
    public int size() {
        return count;
    }

    @Override
    public void swapRemove(int i) {
        if (objects[count - 1] != null) {
            if (i != count - 1) {
                objects[i] = objects[count - 1];
            }
            this.objects[count - 1] = null;
        } else {
            objects[i] = null;
        }
        this.count--;
    }

}
