package etomica.atom;

import etomica.space.CoordinateFactory;
import etomica.species.Species;

/**
 * Builder of a monoatomic atom group, which comprises just an Atom.
 *
 * @author David Kofke
 */

public class AtomFactoryMono extends AtomFactory {
    
    public AtomFactoryMono(CoordinateFactory coordFactory, AtomTypeLeaf atomType, AtomSequencerFactory seqFactory) {
        super(coordFactory, atomType, seqFactory, AtomTreeNodeLeaf.FACTORY);
    }
    
    public void setSpecies(Species species) {
        atomType.setSpecies(species);
    }

    /**
     * Returns a new leaf atom having no children.
     */
    public Atom makeAtom() {
        return newParentAtom();
    }
    
}//end of AtomFactoryMono