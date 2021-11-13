package etomica.util.collections;

import etomica.util.Debug;

import java.util.AbstractList;
import java.util.Arrays;
import java.util.RandomAccess;

public final class IntArrayList {
    private int[] data;
    private int size;

    public IntArrayList(int initialCapacity) {
        super();
        this.data = new int[initialCapacity];
        this.size = 0;
    }

    private void grow() {
        this.data = Arrays.copyOf(this.data, this.data.length * 2);
//        System.out.println("grow " + data.length);
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

    public void clear() {
        this.size = 0;
    }

//    @Override
//    public Integer get(int i) {
//        return getInt(i);
//    }
//
//    @Override
    public int size() {
        return size;
    }
}
