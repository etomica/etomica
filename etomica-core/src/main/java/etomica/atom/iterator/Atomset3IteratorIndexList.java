/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.atom.iterator;


import etomica.atom.AtomListFromArray;
import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.molecule.IMolecule;
import etomica.molecule.IMoleculeList;

/**
 * Atomset Iterator that iterates over set-of-three atoms
 *
 * @author Tai Tan
 */
public class Atomset3IteratorIndexList implements AtomsetIteratorBasisDependent {

    private final int[][] index;
    private final AtomListFromArray atomset;
    private final IAtom[] atoms;
    private int cursor;
    private IMolecule parentGroup;
    private IAtom target;

    /**
     * Constructs iterator without defining set of atoms.
     */
    public Atomset3IteratorIndexList(int[][] index) {
        atomset = new AtomListFromArray(3);
        atoms = atomset.getArray();
        this.index = index;
    }

    public int basisSize() {
        return 1;
    }

    public boolean haveTarget(IAtom a) {
        if (parentGroup == null) {
            return false;
        }

        for (int i = 0; i < index.length; i++) {   //index.length = number of sets

            if (a == parentGroup.getChildList().get(index[i][0])) {
                return true;
            }
            if (a == parentGroup.getChildList().get(index[i][1])) {
                return true;
            }
            if (a == parentGroup.getChildList().get(index[i][2])) {
                return true;
            }
        }
        return false;
    }

    public void setTarget(IAtom a) {
        target = a;
        unset();
    }

    public void setBasis(IMoleculeList parent) {
        if (parent == null) {
            parentGroup = null;
        } else {
            parentGroup = parent.get(0);
        }
        unset();
    }

    public int size() {
        return index.length;
    }

    /**
     * Returns true if three non-null atoms have set and a call to reset() has
     * been performed, without any subsequent calls to next() or nextPair().
     */
    protected boolean hasNext() {

        if (target != null) {
            for (; cursor < index.length; cursor++) {   //index.length = number of pairs

                if (target == parentGroup.getChildList().get(index[cursor][0])) {
                    break;
                }
                if (target == parentGroup.getChildList().get(index[cursor][1])) {
                    break;
                }
                if (target == parentGroup.getChildList().get(index[cursor][2])) {
                    break;
                }
            }
        }
        return (cursor < index.length);

    }

    /**
     * Sets iterator to a state where hasNext() returns false.
     */
    public void unset() {
        cursor = index.length;
    }

    /**
     * Resets iterator to a state where hasNext is true, if atoms in pair are
     * not null.
     */
    public void reset() {
        if (parentGroup == null) {
            return;
        }
        cursor = 0;
    }

    /**
     * Returns the iterator's pair and unsets iterator.
     */
    public AtomListFromArray nextSet() {
        if (!hasNext())
            return null;
        atoms[0] = parentGroup.getChildList().get(index[cursor][0]);
        atoms[1] = parentGroup.getChildList().get(index[cursor][1]);
        atoms[2] = parentGroup.getChildList().get(index[cursor][2]);

        cursor++;
        return atomset;
    }

    /**
     * Same as nextSet().
     */
    public IAtomList next() {
        return nextSet();
    }

    /**
     * Returns 3, indicating that this is a set of three iterator.
     */
    public final int nBody() {
        return 3;
    }

}

