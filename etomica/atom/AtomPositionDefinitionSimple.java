package etomica.atom;

import etomica.space.Vector;


/**
 * simply returns (leaf) atom's position 
 */
public class AtomPositionDefinitionSimple implements AtomPositionDefinition, java.io.Serializable {

    public Vector position(Atom atom) {
        return atom.coord.position();
    }

}
