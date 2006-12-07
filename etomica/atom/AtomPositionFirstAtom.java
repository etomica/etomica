package etomica.atom;

import etomica.space.Vector;

/**
 * Returns the position of the first child leaf atom.  Recurses to find
 * the first child leaf atom.
 */

public class AtomPositionFirstAtom implements AtomPositionDefinition, java.io.Serializable {

    public Vector position(Atom atom) {
        return atom.getNode().firstLeafAtom().coord.position();
    }

    private static final long serialVersionUID = 1L;
}
