package etomica.atom;

import etomica.Atom;
import etomica.AtomFactory;
import etomica.AtomType;
import etomica.Space;
import etomica.Species;

/**
 * Builder of a monoatomic atom group, which comprises just an Atom.
 *
 * @author David Kofke
 */

public class AtomFactoryMono extends AtomFactory {
    
    /**
     * Constructor with AtomTypeSphere default.
     */
    public AtomFactoryMono(Space space, AtomSequencerFactory seqFactory) {
        this(space, new AtomTypeSphere(), seqFactory);
    }
    
    public AtomFactoryMono(Space space, AtomType atomType, AtomSequencerFactory seqFactory) {
        super(space, atomType, seqFactory, AtomTreeNodeLeaf.FACTORY);
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