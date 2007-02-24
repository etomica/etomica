package etomica.atom;

import etomica.space.IVector;


/**
 * Returns a vector given an atom, thereby defining the position
 * of the atom or atom group.  Example implementations of this interface
 * are based on the center of mass, or on the position of the first
 * leaf atom in the group.
 */

/*
 * History
 * Created on Jan 27, 2005 by kofke
 */
public interface AtomPositionDefinition {

    /**
     * Returns the defined position for the given atom, which 
     * may be an atom group.
     */
    public IVector position(Atom atom);
}
