package etomica.atom;

import etomica.api.IAtomKinetic;
import etomica.api.IAtomTypeLeaf;
import etomica.api.IVectorMutable;
import etomica.space.ISpace;

public class AtomLeafDynamic extends AtomLeaf implements IAtomKinetic {

    public AtomLeafDynamic(ISpace space, IAtomTypeLeaf type) {
        super(space, type);
        velocity = space.makeVector();
    }
    
    public IVectorMutable getVelocity() {
        return velocity;
    }
    
    private static final long serialVersionUID = 1L;
    protected final IVectorMutable velocity;
}
