package etomica.box.storage;

import etomica.atom.IAtom;
import etomica.molecule.IMolecule;

import java.util.Arrays;
import java.util.BitSet;

public class IntStorage implements Storage {
    private static final double SIZE_INCREASE_RATIO = 0.3;
    private int[] data;
//    private DoubleWrapper[] views;
    private final BitSet validBits;
    private int count;

    public IntStorage() {
        this.count = 0;
        this.data = new int[count];
//        this.views = new DoubleWrapper[count];
        this.validBits = new BitSet(count);
    }

    public int get(int i) {
        return this.data[i];
    }

//    public double get(IAtom atom) {
//        return this.get(atom.getLeafIndex());
//    }
//
//    public double get(IMolecule mol) {
//        return this.get(mol.getGlobalIndex());
//    }
//
    public void create(int i, int value) {
        validBits.set(i);
        this.set(i, value);
    }

    public void set(int i, int value) {
        this.data[i] = value;
    }

    private void resize(int newSize) {
        this.data = Arrays.copyOf(data, newSize);
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
        this.validBits.set(i, validBits.get(count - 1));
        this.validBits.clear(count - 1);
        count--;
    }
}
