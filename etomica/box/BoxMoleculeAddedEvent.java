package etomica.box;

import etomica.api.IBox;
import etomica.api.IBoxMoleculeEvent;
import etomica.api.IMolecule;


/**
 * Event that conveys that an Atom has been added to a Box.
 */
public class BoxMoleculeAddedEvent extends BoxMoleculeEvent implements IBoxMoleculeEvent {

    public BoxMoleculeAddedEvent(IBox box, IMolecule mole) {
        super(box, mole);
    }

    private static final long serialVersionUID = 1L;
}
