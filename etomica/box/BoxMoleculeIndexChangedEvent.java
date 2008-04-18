package etomica.box;

import etomica.api.IBox;
import etomica.api.IMolecule;


/**
 * Event that conveys that an Atom's global index in a box has changed.
 */
public class BoxMoleculeIndexChangedEvent extends BoxAtomEvent {

    public BoxMoleculeIndexChangedEvent(IBox box, IMolecule atom, int oldIndex) {
        super(box, atom);
        this.oldIndex = oldIndex;
    }
    
    public int getOldIndex() {
        return oldIndex;
    }
    
    private final int oldIndex;
    private static final long serialVersionUID = 1L;
}
