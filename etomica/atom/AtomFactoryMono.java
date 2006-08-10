package etomica.atom;

import etomica.space.CoordinateFactory;
import etomica.species.Species;

/**
 * Builder of a monoatomic atom group, which comprises just an Atom.
 *
 * @author David Kofke
 */

public class AtomFactoryMono extends AtomFactory {
    
    public AtomFactoryMono(CoordinateFactory coordFactory, AtomTypeLeaf atomType) {
        super(atomType, AtomTreeNodeLeaf.FACTORY);
        this.coordFactory = coordFactory;
    }

    /**
     * Returns the CoordinateFactory that gives coordinates to the
     * atom (or the root atom, if this makes an atom group) made by this
     * AtomFactory
     */
    public CoordinateFactory getCoordinateFactory() {
        return coordFactory;
    }

    /**
     * Returns a new leaf atom having no children.
     */
    public Atom makeAtom() {
        return new AtomLeaf(coordFactory.makeCoordinate(), atomType, nodeFactory);
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
    
    protected final CoordinateFactory coordFactory;
    
}//end of AtomFactoryMono