package etomica.atom;

import etomica.Atom;
import etomica.Space;

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
        super(space, seqFactory, AtomTreeNodeLeaf.FACTORY);
        setType(new AtomTypeSphere(this));//default
    }
    
    public boolean isGroupFactory() {return false;}
    
    public void setType(AtomType t) {atomType = t;}
    public AtomType type() {return atomType;}
    
    /**
     * Returns a new leaf atom having no children.
     */
    public Atom makeAtom() {
        return newParentAtom();
    }
    
}//end of AtomFactoryMono