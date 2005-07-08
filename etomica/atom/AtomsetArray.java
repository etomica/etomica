package etomica.atom;

import etomica.Atom;
import etomica.AtomSet;

/**
 * AtomSet formed by wrapping an Atom array.  Size of array
 * cannot be changed after construction.
 */

/*
 * History Created on Feb 18, 2005 by kofke
 */
public class AtomsetArray implements AtomSet, java.io.Serializable {

    /**
     * Wraps a new atom array of the given length.
     */
    public AtomsetArray(int nAtoms) {
        atoms = new Atom[nAtoms];
    }

    /**
     * Makes a new instance holding the atoms in the given atom set. Makes
     * zero-body AtomSet if argument is null.
     */
    public AtomsetArray(AtomSet atomSet) {
        this((atomSet != null) ? atomSet.count() : 0);
        for (int i = 0; i < atoms.length; i++) {
            atoms[i] = atomSet.getAtom(i);
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
        this((Atom[]) atomSet.atoms.clone());

    }

    /**
     * Wraps the given atom array. Subsequent call to getArray will return the
     * array instance given here.
     */
    public AtomsetArray(Atom[] atoms) {
        this.atoms = atoms;
    }

    /**
     * Part of implementation of AtomSet interface.
     */
    public Atom getAtom(int i) {
        return atoms[i];
    }

    /**
     * @return the wrapped array of atoms, which is declared final in the class.
     */
    public Atom[] getArray() {
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
    public void setAtoms(Atom[] newAtoms) {
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
    public void setAtoms(AtomSet atomSet) {
        if (atomSet.count() != atoms.length)
            throw new IllegalArgumentException("Wrong size for atomSet");
        for (int i = 0; i < atoms.length; i++) {
            atoms[i] = atomSet.getAtom(i);
        }
    }

    /**
     * Returns the length of the wrapped array.
     */
    public int count() {
        return atoms.length;
    }

    public boolean equals(Object object) {
        if (!(object instanceof AtomSet)) {
            return false;
        }
        return equals((AtomSet) object);
    }

    /**
     * Returns true if elements of the wrapped array are equal to those in the
     * given atom set, with the comparison performed element-by-element.  Returns
     * false otherwise, or if argument is null, or if atoms.count() != count().
     */
    public boolean equals(AtomSet atoms) {
        if (atoms == null || atoms.count() != count()) {
            return false;
        }
        for (int i = 0; i < this.atoms.length; i++) {
            if (this.atoms[i] != atoms.getAtom(i))
                return false;
        }
        return true;
    }

    private final Atom[] atoms;
}
