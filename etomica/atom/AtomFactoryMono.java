package etomica.atom;

import etomica.Atom;
import etomica.AtomFactory;
import etomica.AtomTreeNodeLeaf;
import etomica.AtomTypeLeaf;
import etomica.Species;
import etomica.space.CoordinateFactory;

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