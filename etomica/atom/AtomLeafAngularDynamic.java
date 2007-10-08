package etomica.atom;

import etomica.space.IVector;
import etomica.space.IOrientation;
import etomica.space.Space;

public class AtomLeafAngularDynamic extends AtomLeafDynamic implements
        IAtomOrientedKinetic {

    private static final long serialVersionUID = 1L;
    public AtomLeafAngularDynamic(Space space, AtomType type) {
        super(space, type);
        iOrientation = space.makeOrientation();
        angularVelocity = space.makeVector();  //XXX wrong! see https://rheneas.eng.buffalo.edu/bugzilla/show_bug.cgi?id=128
    }

    public IVector getAngularVelocity() {
        return angularVelocity;
    }

    public IOrientation getOrientation() {
        return iOrientation;
    }

    protected final IOrientation iOrientation;
    protected final IVector angularVelocity;
}
