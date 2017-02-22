package etomica.atom;

import etomica.api.IAtomType;
import etomica.api.IVectorMutable;
import etomica.space.ISpace;
import etomica.spaceNd.VectorND;

public class AtomOrientedQuaternion extends Atom {

    protected final VectorND quaternion;
    
    public AtomOrientedQuaternion(ISpace space, IAtomType type) {
        super(space, type);
        quaternion = new VectorND(4);
        quaternion.setX(0, 1);
    }

    public IVectorMutable getQuaternion() {
        return quaternion;
    }
}
