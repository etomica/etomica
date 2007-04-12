package etomica.atom;

import etomica.space.Space;

/**
 * Builder of a monoatomic atom group, which comprises just an Atom.
 *
 * @author David Kofke
 */

public class AtomFactoryMono extends AtomFactory {
    
    public AtomFactoryMono(Space space, AtomTypeLeaf atomType) {
        super(atomType);
        this.space = space;
    }

    /**
     * Returns the CoordinateFactory that gives coordinates to the
     * atom (or the root atom, if this makes an atom group) made by this
     * AtomFactory
     */
    public Space getSpace() {
        return space;
    }
    
    /**
     * Returns a new leaf atom having no children.
     */
    public IAtom makeAtom() {
        isMutable = false;
        return new AtomLeaf(space, atomType);
    }
    
    /**
     * Returns 1, indicating that each atom produced by this factory is a single
     * atom with no children.
     */
    public int getNumTreeAtoms() {
        return 1;
    }

    /**
     * Returns 0, becuase this factory makes a leaf atoms, having no children.
     */
    public int getNumChildAtoms() {
        return 0;
    }

    /**
     * Returns 1.
     */
    public int getNumLeafAtoms() {
        return 1;
    }
    
    private static final long serialVersionUID = 1L;
    protected Space space;
}
