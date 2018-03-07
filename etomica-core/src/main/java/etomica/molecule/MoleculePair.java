/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.molecule;


import java.util.*;

/**
 * Data structure that contains two mutable atom instances.
 */
public final class MoleculePair extends AbstractList<IMolecule> implements IMoleculeList, RandomAccess {
    public IMolecule mol0, mol1;

    public MoleculePair() {
    }

    public MoleculePair(IMolecule mol0, IMolecule mol1) {
        this.mol0 = mol0;
        this.mol1 = mol1;
    }

    public boolean equals(Object obj) {
        if (((MoleculePair)obj).mol0 == mol0) {
            return ((MoleculePair)obj).mol1 == mol1;
        }
        return (((MoleculePair)obj).mol0 == mol1 && ((MoleculePair)obj).mol1 == mol0);
    }

    public int hashCode() {
        return mol0.hashCode() + mol1.hashCode();
    }

    public IMolecule get(int i) {
        if(i == 0) return mol0;
        if(i == 1) return mol1;
        throw new IllegalArgumentException();
    }

    @Override
    public IMolecule set(int i, IMolecule iAtom) {
        IMolecule oldAtom;
        switch (i) {
            case 0:
                oldAtom = mol0;
                mol0 = iAtom;
                break;
            case 1:
                oldAtom = mol1;
                mol1 = iAtom;
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
        if (o == mol0) {
            return 0;
        } else if (o == mol1) {
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
        return mol0 == null && mol1 == null;
    }

    @Override
    public boolean contains(Object o) {
        return o == mol0 || o == mol1;
    }

    @Override
    public Iterator<IMolecule> iterator() {
        return new Itr();
    }

    @Override
    public IMolecule[] toArray() {
        return new IMolecule[]{mol0, mol1};
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
        mol0 = null;
        mol1 = null;
    }

    public String toString() {
        return "["+ mol0 +","+ mol1 +"]";
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
                    return mol0;
                case 1:
                    cursor++;
                    return mol1;
                default:
                    throw new NoSuchElementException();

            }
        }
    }

}
