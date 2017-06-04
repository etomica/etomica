package etomica.atom;

import etomica.space.Vector;
import etomica.space.Space;
import etomica.spaceNd.VectorND;

public class AtomOrientedQuaternion extends Atom {

    protected final VectorND quaternion;
    
    public AtomOrientedQuaternion(Space space, IAtomType type) {
        super(space, type);
        quaternion = new VectorND(4);
        quaternion.setX(0, 1);
    }

    public Vector getQuaternion() {
        return quaternion;
    }
}
