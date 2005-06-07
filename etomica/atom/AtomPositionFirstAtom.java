package etomica.atom;

import etomica.Atom;
import etomica.space.Vector;

/**
 * Returns the position of the first child leaf atom.  Recurses to find
 * the first child leaf atom.
 */

public class AtomPositionFirstAtom implements AtomPositionDefinition {

    public Vector position(Atom atom) {
        return atom.node.firstLeafAtom().coord.position();
    }

}