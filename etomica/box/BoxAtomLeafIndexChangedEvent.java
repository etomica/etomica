package etomica.box;

import etomica.api.IAtom;
import etomica.api.IBox;
import etomica.api.IBoxAtomLeafIndexChangedEvent;


/**
 * Event that conveys that an Atom's global index in a box has changed.
 */
public class BoxAtomLeafIndexChangedEvent extends BoxAtomEvent implements IBoxAtomLeafIndexChangedEvent {

    public BoxAtomLeafIndexChangedEvent(IBox box, IAtom atom, int oldIndex) {
        super(box, atom);
        this.oldIndex = oldIndex;
    }
    
    /* (non-Javadoc)
     * @see etomica.box.IBoxAtomLeafIndexChangedEvent#getOldIndex()
     */
    public int getOldIndex() {
        return oldIndex;
    }
    
    private final int oldIndex;
    private static final long serialVersionUID = 1L;
}
