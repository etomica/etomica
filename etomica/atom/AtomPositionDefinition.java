package etomica.atom;

import etomica.api.IAtom;
import etomica.api.IVector;


/**
 * Returns a vector given an atom, thereby defining the position
 * of the atom or atom group.  Example implementations of this interface
 * are based on the center of mass, or on the position of the first
 * leaf atom in the group.
 */
public interface AtomPositionDefinition {

    /**
     * Returns the defined position for the given atom, which 
     * may be an atom group.
     */
    public IVector position(IAtom atom);
}
