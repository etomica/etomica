package etomica.box;

import etomica.api.IBox;
import etomica.api.IMolecule;


/**
 * Event that conveys that an Atom's global index in a box has changed.
 */
public class BoxMoleculeIndexChangedEvent extends BoxMoleculeEvent {

    public BoxMoleculeIndexChangedEvent(IBox box, IMolecule mole, int oldIndex) {
        super(box, mole);
        this.oldIndex = oldIndex;
    }
    
    public int getOldIndex() {
        return oldIndex;
    }
    
    private final int oldIndex;
    private static final long serialVersionUID = 1L;
}
