package etomica.atom;

import etomica.space.ICoordinateKinetic;
import etomica.space.IVector;
import etomica.space.Space;

public class AtomLeafDynamic extends AtomLeaf implements ICoordinateKinetic {

    public AtomLeafDynamic(Space space, AtomType type) {
        super(space, type);
        velocity = space.makeVector();
    }
    
    public IVector getVelocity() {
        return velocity;
    }
    
    private static final long serialVersionUID = 1L;
    protected final IVector velocity;
}
