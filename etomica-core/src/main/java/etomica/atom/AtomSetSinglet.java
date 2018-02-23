/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.atom;


import java.util.*;

/**
 * Data structure that contains a single mutable atom instance.
 */
public final class AtomSetSinglet extends AbstractList<IAtom> implements IAtomList {

    public IAtom atom;

    public AtomSetSinglet() {
    }

    public AtomSetSinglet(IAtom atom) {
        this.atom = atom;
    }

    public IAtom get(int i) {
        if (i == 0) return atom;
        throw new IllegalArgumentException();
    }

    @Override
    public IAtom set(int i, IAtom iAtom) {
        if (i == 0) {
            atom = iAtom;
        }
        throw new IndexOutOfBoundsException();
    }

    @Override
    public void add(int i, IAtom iAtom) {
        throw new UnsupportedOperationException();
    }

    @Override
    public IAtom remove(int i) {
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
    public Iterator<IAtom> iterator() {
        return new Itr();
    }

    @Override
    public IAtom[] toArray() {
        return new IAtom[]{atom};
    }

    @Override
    public boolean add(IAtom iAtom) {
        throw new UnsupportedOperationException();
    }

    @Override
    public boolean remove(Object o) {
        throw new UnsupportedOperationException();
    }

    @Override
    public boolean addAll(Collection<? extends IAtom> collection) {
        throw new UnsupportedOperationException();
    }

    @Override
    public boolean addAll(int i, Collection<? extends IAtom> collection) {
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

    private class Itr implements Iterator<IAtom> {
        private boolean didNext = false;

        @Override
        public boolean hasNext() {
            return !didNext;
        }

        @Override
        public IAtom next() {
            if (!didNext) {
                didNext = true;
                return atom;
            }
            throw new NoSuchElementException();
        }
    }
}
