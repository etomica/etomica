package etomica.atom;

import etomica.space.ICoordinateAngularKinetic;
import etomica.space.IVector;
import etomica.space.Orientation;
import etomica.space.Space;

public class AtomLeafAngularDynamic extends AtomLeafDynamic implements
        ICoordinateAngularKinetic {

    private static final long serialVersionUID = 1L;
    public AtomLeafAngularDynamic(Space space, AtomType type) {
        super(space, type);
        orientation = space.makeOrientation();
        angularVelocity = space.makeVector();  //XXX wrong! see https://rheneas.eng.buffalo.edu/bugzilla/show_bug.cgi?id=128
    }

    public IVector getAngularVelocity() {
        return angularVelocity;
    }

    public Orientation getOrientation() {
        return orientation;
    }

    protected final Orientation orientation;
    protected final IVector angularVelocity;
}
