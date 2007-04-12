package etomica.atom;

import etomica.space.Space;

public class AtomFactoryMonoDynamic extends AtomFactoryMono {

    private static final long serialVersionUID = 1L;

    public AtomFactoryMonoDynamic(Space space, AtomTypeLeaf atomType) {
        super(space, atomType);
    }

    /**
     * Returns a new leaf atom having no children.
     */
    public IAtom makeAtom() {
        isMutable = false;
        return new AtomLeafDynamic(space, atomType);
    }
    
}
