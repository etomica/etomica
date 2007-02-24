package etomica.atom;

import etomica.space.Vector;


/**
 * simply returns (leaf) atom's position 
 */
public class AtomPositionDefinitionSimple implements AtomPositionDefinition, java.io.Serializable {

    public Vector position(Atom atom) {
        return ((AtomLeaf)atom).getCoord().getPosition();
    }

    private static final long serialVersionUID = 1L;
}
