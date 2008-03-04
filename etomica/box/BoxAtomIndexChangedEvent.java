package etomica.box;

import etomica.api.IAtom;
import etomica.api.IBox;


/**
 * Event that conveys that an Atom's global index in a box has changed.
 */
public class BoxAtomIndexChangedEvent extends BoxAtomEvent {

    public BoxAtomIndexChangedEvent(IBox box, IAtom atom, int oldIndex) {
        super(box, atom);
        this.oldIndex = oldIndex;
    }
    
    public int getOldIndex() {
        return oldIndex;
    }
    
    private final int oldIndex;
    private static final long serialVersionUID = 1L;
}
