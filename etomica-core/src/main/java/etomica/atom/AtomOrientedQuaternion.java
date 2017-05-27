package etomica.atom;

import etomica.api.IAtomType;
import etomica.api.IVector;
import etomica.space.Space;
import etomica.spaceNd.VectorND;

public class AtomOrientedQuaternion extends Atom {

    protected final VectorND quaternion;
    
    public AtomOrientedQuaternion(Space space, IAtomType type) {
        super(space, type);
        quaternion = new VectorND(4);
        quaternion.setX(0, 1);
    }

    public IVector getQuaternion() {
        return quaternion;
    }
}
