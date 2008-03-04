package etomica.box;

import etomica.api.IAtom;
import etomica.api.IBox;

/**
 * Event that conveys some happening with respect to an Atom in a Box.
 */
public class BoxAtomEvent extends BoxEvent {
    
    public BoxAtomEvent(IBox box, IAtom atom) {
        super(box);
        this.atom = atom;
    }

    public IAtom getAtom() {
        return atom;
    }
    
    private final IAtom atom;
    private static final long serialVersionUID = 1L;
}
