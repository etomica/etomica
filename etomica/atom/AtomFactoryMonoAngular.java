package etomica.atom;

import etomica.space.Space;

public class AtomFactoryMonoAngular extends AtomFactoryMono {

    private static final long serialVersionUID = 1L;

    public AtomFactoryMonoAngular(Space space, AtomTypeOrientedSphere atomType) {
        super(space, atomType);
    }

    /**
     * Returns a new leaf atom having no children.
     */
    public IAtom makeAtom() {
        isMutable = false;
        return new AtomLeafAngular(space, atomType);
    }
    
}
