package etomica.box;

import etomica.api.IBox;
import etomica.api.IBoxMoleculeRemovedEvent;
import etomica.api.IMolecule;


/**
 * Event that conveys that an Atom has been removed from a Box.
 */
public class BoxMoleculeRemovedEvent extends BoxMoleculeEvent implements IBoxMoleculeRemovedEvent {

    public BoxMoleculeRemovedEvent(IBox box, IMolecule mole) {
        super(box, mole);
    }

    private static final long serialVersionUID = 1L;
}
