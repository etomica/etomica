package etomica.atom;

import etomica.api.IAtom;
import etomica.api.IAtomPositioned;
import etomica.api.IVector;


/**
 * simply returns (leaf) atom's position 
 */
public class AtomPositionDefinitionSimple implements AtomPositionDefinition, java.io.Serializable {

    public IVector position(IAtom atom) {
        return ((IAtomPositioned)atom).getPosition();
    }

    private static final long serialVersionUID = 1L;
}
