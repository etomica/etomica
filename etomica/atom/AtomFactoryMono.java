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
    
    public void setSpecies(Species species) {
        atomType.setSpecies(species);
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
    
    protected final CoordinateFactory coordFactory;
    
}//end of AtomFactoryMono