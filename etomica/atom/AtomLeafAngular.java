package etomica.atom;

import etomica.space.IOrientation;
import etomica.space.Space;

public class AtomLeafAngular extends AtomLeaf implements
        IAtomOriented {

    public AtomLeafAngular(Space space, AtomTypeLeaf type) {
        super(space, type);
        iOrientation = space.makeOrientation();
    }

    public IOrientation getOrientation() {
        return iOrientation;
    }

    private static final long serialVersionUID = 1L;
    protected final IOrientation iOrientation;
}
