/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.atom;


/**
 * AtomSet formed by wrapping an AtomArrayList.  ArrayList can be 
 */
public class MoleculeListWrapper implements IMoleculeList, java.io.Serializable {

    /**
     * Wraps a new atom array of the given length.
     */
    public MoleculeListWrapper() {
        atoms = new MoleculeArrayList();
    }

    /**
     * Makes a new instance holding the atoms in the given atom set. Makes
     * zero-body AtomSet if argument is null.
     */
    public MoleculeListWrapper(IMoleculeList atomSet) {
        this();
        atoms.ensureCapacity(atomSet.getMoleculeCount());
        for (int i = 0; i < atoms.getMoleculeCount(); i++) {
            atoms.add(atomSet.getMolecule(i));
        }
    }

    /**
     * Wraps the given atom array. Subsequent call to getArray will return the
     * array instance given here.
     */
    public MoleculeListWrapper(MoleculeArrayList atoms) {
        this.atoms = atoms;
    }

    /**
     * Part of implementation of AtomSet interface.
     */
    public IMolecule getMolecule(int i) {
        return atoms.getMolecule(i);
    }

    /**
     * @return the wrapped array of atoms, which is declared final in the class.
     */
    public MoleculeArrayList getArrayList() {
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
    public void setAtoms(IMoleculeList atomSet) {
        atoms.clear();
        atoms.ensureCapacity(atomSet.getMoleculeCount());
        for (int i = 0; i < atomSet.getMoleculeCount(); i++) {
            atoms.add(atomSet.getMolecule(i));
        }
    }

    /**
     * Returns the length of the wrapped array.
     */
    public int getMoleculeCount() {
        return atoms.getMoleculeCount();
    }

    private static final long serialVersionUID = 1L;
    private final MoleculeArrayList atoms;
}
