package etomica.atom;

import etomica.space.IVector;
import etomica.space.Space;


/**
 * Atom type for a rigid molecule with an orientation and (therefore) a moment
 * of inertia.  The molecule holds the orientation and the type holds the
 * moment.  The type does not actually know how to calculate the moment, so
 * the moment must be initialized by whoever creates this type.
 */
public class AtomTypeMoleculeOriented extends AtomTypeMolecule {
    
    public AtomTypeMoleculeOriented(Space space) {
        super(new AtomPositionCOM(space));
        moment = space.makeVector();
    }

    /**
     * Returns the principle components of the moment of inertia of the
     * molecule within the body-fixed frame.
     */
    public IVector momentOfInertia() {
        return moment;
    }

    private static final long serialVersionUID = 1L;
    protected final IVector moment;
}