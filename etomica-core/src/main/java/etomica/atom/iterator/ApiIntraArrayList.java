/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.atom.iterator;

import etomica.atom.IAtomList;
import etomica.atom.AtomArrayList;
import etomica.atom.AtomPair;
import etomica.util.Debug;

/**
 * Returns all pairs formed from a single list of atoms.
 */
public class ApiIntraArrayList implements AtomLeafsetIterator, java.io.Serializable {

    /**
     * Construct iterator with an empty list. No iterates will be given until a
     * non-empty list is specified via setList.
     */
    public ApiIntraArrayList() {
        this(new AtomArrayList());
    }

    /**
     * Constructs iterator to return pairs from the given list. Requires reset()
     * before first use.
     * 
     * @param list
     */
    public ApiIntraArrayList(IAtomList list) {
        setList(list);
    }

    /**
     * Sets iterator in condition to begin iteration.
     */
    public void reset() {
        if (list.getAtomCount() < 2) {
            outerIndex = 2;
            innerIndex = 2;
            return;
        }
        outerIndex = 0;
        innerIndex = 0;
        atoms.atom0 = list.getAtom(0);
    }

    /**
     * Sets iterator such that next is null.
     */
    public void unset() {
        outerIndex = list.getAtomCount() - 2;
        innerIndex = list.getAtomCount() - 1;
    }

    /**
     * Returns the number of iterates, which is list.size*(list.size-1)/2
     */
    public int size() {
        return list.getAtomCount() * (list.getAtomCount() - 1) / 2;
    }

    /**
     * Returns the next iterate pair. Returns null if hasNext() is false.
     */
    public IAtomList next() {
        if (innerIndex > list.getAtomCount() - 2) {
            if (outerIndex > list.getAtomCount() - 3) {
                return null;
            }
            outerIndex++;
            atoms.atom0 = list.getAtom(outerIndex);
            innerIndex = outerIndex;
        }
        innerIndex++;
        atoms.atom1 = list.getAtom(innerIndex);
        if (Debug.ON && atoms.atom0 == atoms.atom1) {
            throw new RuntimeException("oops");
        }
        return atoms;
    }

    /**
     * Returns 2, indicating that this is a pair iterator
     */
    public int nBody() {
        return 2;
    }

    /**
     * Sets the list that will be used to generate the pairs. Must call reset()
     * before beginning iteration.
     * 
     * @param newList
     *            the new atom list for iteration
     */
    public void setList(IAtomList newList) {
        if (newList.getAtomCount() > 1 && newList.getAtom(0) == newList.getAtom(1)) {
            throw new RuntimeException("oops");
        }
        list = newList;
        unset();
    }

    /**
     * Returns the list used to generate the pairs.
     */
    public IAtomList getList() {
        return list;
    }

    private static final long serialVersionUID = 1L;
    private IAtomList list;
    private int outerIndex, innerIndex;
    private final AtomPair atoms = new AtomPair();
}
