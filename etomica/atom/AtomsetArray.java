package etomica.atom;

import etomica.Atom;
import etomica.AtomSet;


/**
 * AtomSet formed by wrapping an Atom array.
 */

/*
 * History
 * Created on Feb 18, 2005 by kofke
 */
public class AtomsetArray implements AtomSet {

    /**
     * Wraps a new atom array of the given length.
     */
    public AtomsetArray(int nAtoms) {
        atoms = new Atom[nAtoms];
    }
    
    /**
     * Wraps the given atom array.  Subsequent call to getArray
     * will return the array instance given here.
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
     */
    public void setAtoms(Atom[] newAtoms) {
        if(newAtoms.length != atoms.length) throw new IllegalArgumentException("Wrong size array");
        System.arraycopy(newAtoms, 0, this.atoms, 0, atoms.length);
    }

    /**
     * Returns the length of the wrapped array.
     */
    public int count() {
        return atoms.length;
    }

    /**
     * Returns true if elements of the wrapped array are equal
     * to those in the given atom set, with the comparison
     * performed element-by-element.
     */
    public boolean equals(AtomSet atoms) {
        if(atoms.count() != count()) return false;
        for(int i=0; i<this.atoms.length; i++) {
            if(this.atoms[i] != atoms.getAtom(i)) return false;
        }
        return true;
    }

    private final Atom[] atoms;
}
