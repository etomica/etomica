package etomica.atom;

import etomica.Atom;
import etomica.AtomSet;


/**
 * AtomSet formed by wrapping an Atom array
 *
 */

/*
 * History
 * Created on Feb 18, 2005 by kofke
 */
public class AtomsetArray implements AtomSet {

    /**
     * 
     */
    public AtomsetArray(int nAtoms) {
        atoms = new Atom[nAtoms];
    }

    /* (non-Javadoc)
     * @see etomica.AtomSet#getAtom(int)
     */
    public Atom getAtom(int i) {
        return atoms[i];
    }
    
    public void setAtoms(Atom[] newAtoms) {
        if(newAtoms.length != atoms.length) throw new IllegalArgumentException("Wrong size array");
        System.arraycopy(newAtoms, 0, this.atoms, 0, atoms.length);
    }

    /* (non-Javadoc)
     * @see etomica.AtomSet#count()
     */
    public int count() {
        return atoms.length;
    }

    /* (non-Javadoc)
     * @see etomica.AtomSet#equals(etomica.AtomSet)
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
