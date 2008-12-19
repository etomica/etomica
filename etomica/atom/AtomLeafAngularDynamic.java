package etomica.atom;

import etomica.api.IAtomTypeLeaf;
import etomica.api.IVectorMutable;
import etomica.space.IOrientation;
import etomica.space.ISpace;

public class AtomLeafAngularDynamic extends AtomLeafDynamic implements
        IAtomOrientedKinetic {

    private static final long serialVersionUID = 1L;
    public AtomLeafAngularDynamic(ISpace space, IAtomTypeLeaf type) {
        super(space, type);
        iOrientation = space.makeOrientation();
        angularVelocity = space.makeVector();  //XXX wrong! see https://rheneas.eng.buffalo.edu/bugzilla/show_bug.cgi?id=128
    }

    public IVectorMutable getAngularVelocity() {
        return angularVelocity;
    }

    public IOrientation getOrientation() {
        return iOrientation;
    }

    protected final IOrientation iOrientation;
    protected final IVectorMutable angularVelocity;
}
