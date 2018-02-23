/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.atom;


/**
 * AtomSet formed by wrapping an Atom array.  Size of array
 * cannot be changed after construction.
 */
public class AtomsetArray implements IAtomList, java.io.Serializable {

    /**
     * Wraps a new atom array of the given length.
     */
    public AtomsetArray(int nAtoms) {
        atoms = new IAtom[nAtoms];
    }

    /**
     * Makes a new instance holding the atoms in the given atom set. Makes
     * zero-body AtomSet if argument is null.
     */
    public AtomsetArray(IAtomList atomSet) {
        this((atomSet != null) ? atomSet.size() : 0);
        for (int i = 0; i < atoms.length; i++) {
            atoms[i] = atomSet.get(i);
        }
    }

    /**
     * Copy constructor, wrapping a new array unique to this instance but
     * holding the same atoms as the given atom set.
     * 
     * @throws NullPointerException
     *             if argument is null
     */
    public AtomsetArray(AtomsetArray atomSet) {
        this(atomSet.atoms.clone());

    }

    /**
     * Wraps the given atom array. Subsequent call to getArray will return the
     * array instance given here.
     */
    public AtomsetArray(IAtom[] atoms) {
        this.atoms = atoms;
    }

    /**
     * Part of implementation of AtomSet interface.
     */
    public IAtom get(int i) {
        return atoms[i];
    }

    /**
     * @return the wrapped array of atoms, which is declared final in the class.
     */
    public IAtom[] getArray() {
        return atoms;
    }

    /**
     * Copies the atoms in the given array to the wrapped, final array of atoms.
     * 
     * @throws IllegalArgumentException
     *             if length of array is not equal to count field of this
     *             instance.
     * @throws NullPointerException
     *             if argument is null.
     */
    public void setAtoms(IAtom[] newAtoms) {
        if (newAtoms.length != atoms.length)
            throw new IllegalArgumentException("Wrong size array; should be "
                    + atoms.length + " but is " + newAtoms.length);
        System.arraycopy(newAtoms, 0, this.atoms, 0, atoms.length);
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
        if (atomSet.size() != atoms.length)
            throw new IllegalArgumentException("Wrong size for atomSet");
        for (int i = 0; i < atoms.length; i++) {
            atoms[i] = atomSet.get(i);
        }
    }

    /**
     * Returns the length of the wrapped array.
     */
    public int size() {
        return atoms.length;
    }
    
    public String toString() {
        String str = "[";
        if (atoms.length > 0) {
            str += atoms[0].toString();
        }
        for (int i=1; i<atoms.length; i++) {
            str += ", "+atoms[i].toString();
        }
        str += "]";
        return str;
    }

    private static final long serialVersionUID = 1L;
    private final IAtom[] atoms;
}
