package etomica.atom;

import etomica.api.IBox;
import etomica.api.IMolecule;

/**
 * Interface for objects when return atoms (meeting some specification)
 * from a box.
 */
public interface MoleculeSource {
    
    /**
     * sets the Box the source should pull Atoms from.
     * Box should not be null
     */
    public void setBox(IBox p);

    /**
     * Returns an atom.  Will return null if there are no appropriate atoms in 
     * the given box.
     */
    public IMolecule getMolecule();
}
