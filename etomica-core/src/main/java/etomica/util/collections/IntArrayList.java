package etomica.util.collections;

import etomica.util.Debug;

import java.util.Arrays;

public final class IntArrayList {
    private int[] data;
    private int size;

    public IntArrayList(int initialCapacity) {
        super();
        this.data = new int[initialCapacity];
        this.size = 0;
    }

    public IntArrayList(int[] elements) {
        super();
        this.data = elements;
        this.size = elements.length;
    }

    private void grow() {
        this.data = Arrays.copyOf(this.data, this.data.length * 2);
    }

    public IntArrayList() {
        this(32);
    }

    public void ensureCapacity(int capacity) {
        if (this.data.length < capacity) {
            this.data = Arrays.copyOf(this.data, capacity);
        }
    }

    public int getInt(int i) {
        if (Debug.ON && i >= size) {
            throw new IndexOutOfBoundsException();
        }

        return data[i];
    }

    public void add(int x) {
        if (size == data.length) {
            this.grow();
        }

        data[size] = x;
        size++;
    }

    public void set(int i, int x) {
        if (Debug.ON && i >= size) throw new IndexOutOfBoundsException();
        data[i] = x;
    }

    public void clear() {
        this.size = 0;
    }

    public int size() {
        return size;
    }
}
