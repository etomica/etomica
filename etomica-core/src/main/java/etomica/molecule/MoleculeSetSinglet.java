/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.molecule;


import java.util.*;

/**
 * Data structure that contains a single mutable atom instance.
 */
public final class MoleculeSetSinglet extends AbstractList<IMolecule> implements IMoleculeList {

    public IMolecule atom;

    public MoleculeSetSinglet() {
    }

    public MoleculeSetSinglet(IMolecule molecule) {
        this.atom = molecule;
    }

    public IMolecule get(int i) {
        if (i == 0) return atom;
        throw new IllegalArgumentException();
    }

    @Override
    public IMolecule set(int i, IMolecule iMolecule) {
        if (i == 0) {
            atom = iMolecule;
        }
        throw new IndexOutOfBoundsException();
    }

    @Override
    public void add(int i, IMolecule iMolecule) {
        throw new UnsupportedOperationException();
    }

    @Override
    public IMolecule remove(int i) {
        throw new UnsupportedOperationException();
    }

    @Override
    public int indexOf(Object o) {
        return o == atom ? 0 : -1;
    }

    @Override
    public int lastIndexOf(Object o) {
        return indexOf(o);
    }

    public int size() {
        return 1;
    }

    @Override
    public boolean isEmpty() {
        return atom == null;
    }

    @Override
    public boolean contains(Object o) {
        return o == atom;
    }

    @Override
    public Iterator<IMolecule> iterator() {
        return new Itr();
    }

    @Override
    public IMolecule[] toArray() {
        return new IMolecule[]{atom};
    }

    @Override
    public boolean add(IMolecule iMolecule) {
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
        atom = null;
    }

    public String toString() {
        return "[" + atom + "]";
    }

    private class Itr implements Iterator<IMolecule> {
        private boolean didNext = false;

        @Override
        public boolean hasNext() {
            return !didNext;
        }

        @Override
        public IMolecule next() {
            if (!didNext) {
                didNext = true;
                return atom;
            }
            throw new NoSuchElementException();
        }
    }
}
