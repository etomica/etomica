/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.atom;


/**
 * AtomSet formed by wrapping an AtomArrayList.  ArrayList can be 
 */
public class AtomListWrapper implements IAtomList, java.io.Serializable {

    /**
     * Wraps a new atom array of the given length.
     */
    public AtomListWrapper() {
        atoms = new AtomArrayList();
    }

    /**
     * Makes a new instance holding the atoms in the given atom set. Makes
     * zero-body AtomSet if argument is null.
     */
    public AtomListWrapper(IAtomList atomSet) {
        this();
        atoms.ensureCapacity(atomSet.getAtomCount());
        for (int i = 0; i < atoms.getAtomCount(); i++) {
            atoms.add(atomSet.getAtom(i));
        }
    }

    /**
     * Wraps the given atom array. Subsequent call to getArray will return the
     * array instance given here.
     */
    public AtomListWrapper(AtomArrayList atoms) {
        this.atoms = atoms;
    }

    /**
     * Part of implementation of AtomSet interface.
     */
    public IAtom getAtom(int i) {
        return atoms.getAtom(i);
    }

    /**
     * @return the wrapped array of atoms, which is declared final in the class.
     */
    public AtomArrayList getArrayList() {
        return atoms;
    }

    /**
     * Copies the atoms in the given atom set to the wrapped array of atoms.
     * 
     * @throws IllegalArgumentException
     *             if length of array is not equal to count field of this
     *             instance.
     * @throws NullPointerException
     *             if argument is null
     */
    public void setAtoms(IAtomList atomSet) {
        atoms.clear();
        atoms.ensureCapacity(atomSet.getAtomCount());
        for (int i = 0; i < atomSet.getAtomCount(); i++) {
            atoms.add(atomSet.getAtom(i));
        }
    }

    /**
     * Returns the length of the wrapped array.
     */
    public int getAtomCount() {
        return atoms.getAtomCount();
    }

    private static final long serialVersionUID = 1L;
    private final AtomArrayList atoms;
}
