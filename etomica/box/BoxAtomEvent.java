package etomica.box;

import etomica.api.IAtomLeaf;
import etomica.api.IBox;
import etomica.api.IBoxAtomEvent;

/**
 * Event that conveys some happening with respect to an Atom in a Box.
 */
public class BoxAtomEvent extends BoxEvent implements IBoxAtomEvent {
    
    public BoxAtomEvent(IBox box, IAtomLeaf atom) {
        super(box);
        this.atom = atom;
    }

    /* (non-Javadoc)
     * @see etomica.box.IBoxAtomEvent#getAtom()
     */
    public IAtomLeaf getAtom() {
        return atom;
    }
    
    private final IAtomLeaf atom;
    private static final long serialVersionUID = 1L;
}
