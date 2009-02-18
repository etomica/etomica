package etomica.box;

import etomica.api.IBox;
import etomica.api.IBoxMoleculeAddedEvent;
import etomica.api.IMolecule;


/**
 * Event that conveys that an Atom has been added to a Box.
 */
public class BoxMoleculeAddedEvent extends BoxMoleculeEvent implements IBoxMoleculeAddedEvent {

    public BoxMoleculeAddedEvent(IBox box, IMolecule mole) {
        super(box, mole);
    }

    private static final long serialVersionUID = 1L;
}
