package etomica.atom;

import etomica.api.IAtomTypeLeaf;
import etomica.api.IVector;
import etomica.space.ISpace;

public class AtomLeafDynamic extends AtomLeaf implements IAtomKinetic {

    public AtomLeafDynamic(ISpace space, IAtomTypeLeaf type) {
        super(space, type);
        velocity = space.makeVector();
    }
    
    public IVector getVelocity() {
        return velocity;
    }
    
    private static final long serialVersionUID = 1L;
    protected final IVector velocity;
}
