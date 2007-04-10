package etomica.atom;

import etomica.space.ICoordinateAngular;
import etomica.space.Orientation;
import etomica.space.Space;

public class AtomLeafAngular extends AtomLeaf implements
        ICoordinateAngular {

    public AtomLeafAngular(Space space, AtomType type) {
        super(space, type);
        orientation = space.makeOrientation();
    }

    public Orientation getOrientation() {
        return orientation;
    }

    private static final long serialVersionUID = 1L;
    protected final Orientation orientation;
}
