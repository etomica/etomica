/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial;

import etomica.molecule.IMolecule;
import etomica.molecule.IMoleculeList;

import java.util.Collection;
import java.util.Iterator;
import java.util.List;
import java.util.ListIterator;

/**
 * Class that holds a set of atom pairs.  Takes a list of atoms in its
 * constructor, and forms an instance of AtomPair for each pair formed in the
 * list.  Each AtomPair instance can be accessed via the getSPair method. 
 */
public class AtomPairSet implements java.io.Serializable {

    /**
     * Constructor for AtomPairSet.
     * @param list The list of atoms for which the set of pairs is formed.
     */
    public AtomPairSet(IMoleculeList list) {
        aPairs = new MoleculePair[list.size()-1][];
        setAtoms(list);
    }
    
    /**
     * Returns atom pair for ith and jth atoms in set.
     */
    public MoleculePair getAPair(int i, int j) {
        if(i==j) throw new IllegalArgumentException("Error: asking for pair formed with both atoms the same");
        return i<j ? aPairs[i][j-i-1] : aPairs[j][i-j-1];
    }

    private void setAtoms(IMoleculeList list) {
        int N = list.size();
        for(int i=0; i<N-1; i++) {
            aPairs[i] = new MoleculePair[N-1-i];
            for(int j=0; j<N-1-i; j++) {
                IMolecule mol0, mol1;
                mol0 = list.get(i);
                mol1 = list.get(i+j+1);
                if(mol0.getType().getIndex()>mol1.getType().getIndex()){//only for mixtures
                    mol0 = mol1;
                    mol1 = list.get(i);
                }
                aPairs[i][j] = new MoleculePair(mol0, mol1);
            }
        }
    }
    
    private final MoleculePair[][] aPairs;

    public static class MoleculePair implements IMoleculeList {
        protected final IMolecule mol0, mol1;

        public MoleculePair(IMolecule mol0, IMolecule mol1) {
            this.mol0 = mol0;
            this.mol1 = mol1;
        }

        @Override
        public IMolecule get(int i) {
            return i==0 ? mol0 : mol1;
        }

        @Override
        public IMolecule set(int index, IMolecule element) {
            throw new RuntimeException("nope");
        }

        @Override
        public void add(int index, IMolecule element) {
            throw new RuntimeException("nope");
        }

        @Override
        public IMolecule remove(int index) {
            throw new RuntimeException("nope");
        }

        @Override
        public int indexOf(Object o) {
            throw new RuntimeException("nope");
        }

        @Override
        public int lastIndexOf(Object o) {
            throw new RuntimeException("nope");
        }

        @Override
        public ListIterator<IMolecule> listIterator() {
            throw new RuntimeException("nope");
        }

        @Override
        public ListIterator<IMolecule> listIterator(int index) {
            throw new RuntimeException("nope");
        }

        @Override
        public List<IMolecule> subList(int fromIndex, int toIndex) {
            throw new RuntimeException("nope");
        }

        @Override
        public int size() {
            return 2;
        }

        @Override
        public boolean isEmpty() {
            return false;
        }

        @Override
        public boolean contains(Object o) {
            throw new RuntimeException("nope");
        }

        @Override
        public Iterator<IMolecule> iterator() {
            throw new RuntimeException("nope");
        }

        @Override
        public Object[] toArray() {
            throw new RuntimeException("nope");
        }

        @Override
        public <T> T[] toArray(T[] a) {
            throw new RuntimeException("nope");
        }

        @Override
        public boolean add(IMolecule iMolecule) {
            throw new RuntimeException("nope");
        }

        @Override
        public boolean remove(Object o) {
            throw new RuntimeException("nope");
        }

        @Override
        public boolean containsAll(Collection<?> c) {
            throw new RuntimeException("nope");
        }

        @Override
        public boolean addAll(Collection<? extends IMolecule> c) {
            throw new RuntimeException("nope");
        }

        @Override
        public boolean addAll(int index, Collection<? extends IMolecule> c) {
            throw new RuntimeException("nope");
        }

        @Override
        public boolean removeAll(Collection<?> c) {
            throw new RuntimeException("nope");
        }

        @Override
        public boolean retainAll(Collection<?> c) {
            throw new RuntimeException("nope");
        }

        @Override
        public void clear() {
            throw new RuntimeException("nope");
        }
    }
}
