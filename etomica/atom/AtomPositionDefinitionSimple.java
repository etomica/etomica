package etomica.atom;

import etomica.Atom;
import etomica.space.Vector;


/**
 * simply returns (leaf) atom's position 
 */
public class AtomPositionDefinitionSimple implements AtomPositionDefinition {

    public Vector position(Atom atom) {
        return atom.coord.position();
    }

}
