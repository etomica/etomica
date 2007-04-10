package etomica.atom;

import etomica.space.IVector;


/**
 * simply returns (leaf) atom's position 
 */
public class AtomPositionDefinitionSimple implements AtomPositionDefinition, java.io.Serializable {

    public IVector position(Atom atom) {
        return ((AtomLeaf)atom).getPosition();
    }

    private static final long serialVersionUID = 1L;
}
