package etomica.box;

import etomica.api.IBox;
import etomica.api.IBoxMoleculeEvent;
import etomica.api.IMolecule;


/**
 * Event that conveys that an Atom has been removed from a Box.
 */
public class BoxMoleculeRemovedEvent extends BoxMoleculeEvent implements IBoxMoleculeEvent {

    public BoxMoleculeRemovedEvent(IBox box, IMolecule mole) {
        super(box, mole);
    }

    private static final long serialVersionUID = 1L;
}
