package etomica.atom;


/**
 * AtomSet formed by wrapping an AtomArrayList.  ArrayList can be 
 */

/*
 * History Created on Feb 18, 2005 by kofke
 */
public class AtomsetArrayList implements AtomSet, java.io.Serializable {

    /**
     * Wraps a new atom array of the given length.
     */
    public AtomsetArrayList() {
        atoms = new AtomArrayList();
    }

    /**
     * Makes a new instance holding the atoms in the given atom set. Makes
     * zero-body AtomSet if argument is null.
     */
    public AtomsetArrayList(AtomSet atomSet) {
        this();
        atoms.ensureCapacity(atomSet.count());
        for (int i = 0; i < atoms.size(); i++) {
            atoms.add(atomSet.getAtom(i));
        }
    }

    /**
     * Wraps the given atom array. Subsequent call to getArray will return the
     * array instance given here.
     */
    public AtomsetArrayList(AtomArrayList atoms) {
        this.atoms = atoms;
    }

    /**
     * Part of implementation of AtomSet interface.
     */
    public Atom getAtom(int i) {
        return atoms.get(i);
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
    public void setAtoms(AtomSet atomSet) {
        atoms.clear();
        atoms.ensureCapacity(atomSet.count());
        for (int i = 0; i < atomSet.count(); i++) {
            atoms.add(atomSet.getAtom(i));
        }
    }

    /**
     * Returns the length of the wrapped array.
     */
    public int count() {
        return atoms.size();
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
    public boolean equals(AtomSet candidateAtoms) {
        if (candidateAtoms == null || candidateAtoms.count() != count()) {
            return false;
        }
        for (int i = 0; i < atoms.size(); i++) {
            if (atoms.get(i) != candidateAtoms.getAtom(i))
                return false;
        }
        return true;
    }

    private final AtomArrayList atoms;
}
