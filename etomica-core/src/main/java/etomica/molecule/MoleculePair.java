/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.molecule;


import java.util.*;

/**
 * Data structure that contains two mutable atom instances.
 */
public final class MoleculePair extends AbstractList<IMolecule> implements IMoleculeList, RandomAccess {
    public IMolecule atom0, atom1;

    public MoleculePair() {
    }

    public MoleculePair(IMolecule atom0, IMolecule atom1) {
        this.atom0 = atom0;
        this.atom1 = atom1;
    }

    public boolean equals(Object obj) {
        if (((MoleculePair)obj).atom0 == atom0) {
            return ((MoleculePair)obj).atom1 == atom1;
        }
        return (((MoleculePair)obj).atom0 == atom1 && ((MoleculePair)obj).atom1 == atom0);
    }

    public int hashCode() {
        return atom0.hashCode() + atom1.hashCode();
    }

    public IMolecule get(int i) {
        if(i == 0) return atom0;
        if(i == 1) return atom1;
        throw new IllegalArgumentException();
    }

    @Override
    public IMolecule set(int i, IMolecule iAtom) {
        IMolecule oldAtom;
        switch (i) {
            case 0:
                oldAtom = atom0;
                atom0 = iAtom;
                break;
            case 1:
                oldAtom = atom1;
                atom1 = iAtom;
                break;
            default:
                throw new IndexOutOfBoundsException();
        }
        return oldAtom;
    }

    @Override
    public void add(int i, IMolecule iAtom) {
        throw new UnsupportedOperationException();
    }

    @Override
    public IMolecule remove(int i) {
        throw new UnsupportedOperationException();
    }

    @Override
    public int indexOf(Object o) {
        if (o == atom0) {
            return 0;
        } else if (o == atom1) {
            return 1;
        } else {
            return -1;
        }
    }

    public int size() {
        return 2;
    }

    @Override
    public boolean isEmpty() {
        return atom0 == null && atom1 == null;
    }

    @Override
    public boolean contains(Object o) {
        return o == atom0 || o == atom1;
    }

    @Override
    public Iterator<IMolecule> iterator() {
        return new Itr();
    }

    @Override
    public IMolecule[] toArray() {
        return new IMolecule[]{atom0, atom1};
    }


    @Override
    public boolean add(IMolecule iAtom) {
        throw new UnsupportedOperationException();
    }

    @Override
    public boolean remove(Object o) {
        throw new UnsupportedOperationException();
    }

    @Override
    public boolean addAll(Collection<? extends IMolecule> collection) {
        throw new UnsupportedOperationException();
    }

    @Override
    public boolean addAll(int i, Collection<? extends IMolecule> collection) {
        throw new UnsupportedOperationException();
    }

    @Override
    public boolean removeAll(Collection<?> collection) {
        throw new UnsupportedOperationException();
    }

    @Override
    public boolean retainAll(Collection<?> collection) {
        throw new UnsupportedOperationException();
    }

    @Override
    public void clear() {
        atom0 = null;
        atom1 = null;
    }

    public String toString() {
        return "["+atom0+","+atom1+"]";
    }

    private class Itr implements Iterator<IMolecule> {
        private byte cursor = 0;

        @Override
        public boolean hasNext() {
            return cursor < 2;
        }

        @Override
        public IMolecule next() {
            switch (cursor) {
                case 0:
                    cursor++;
                    return atom0;
                case 1:
                    cursor++;
                    return atom1;
                default:
                    throw new NoSuchElementException();

            }
        }
    }

}
