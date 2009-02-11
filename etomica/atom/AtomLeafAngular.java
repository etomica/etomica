package etomica.atom;

import etomica.api.IAtomType;
import etomica.space.IOrientation;
import etomica.space.ISpace;

public class AtomLeafAngular extends AtomLeaf implements
        IAtomOriented {

    public AtomLeafAngular(ISpace space, IAtomType type) {
        super(space, type);
        iOrientation = space.makeOrientation();
    }

    public IOrientation getOrientation() {
        return iOrientation;
    }

    private static final long serialVersionUID = 1L;
    protected final IOrientation iOrientation;
}
