package etomica.atom;

import etomica.space.IVector;

/**
 * Returns the position of the first child leaf atom.  Recurses to find
 * the first child leaf atom.
 */

public class AtomPositionFirstAtom implements AtomPositionDefinition, java.io.Serializable {

    public IVector position(Atom atom) {
        return atom.firstLeafAtom().getCoord().getPosition();
    }

    private static final long serialVersionUID = 1L;
}
