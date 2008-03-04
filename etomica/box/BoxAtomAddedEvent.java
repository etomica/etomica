package etomica.box;

import etomica.api.IAtom;
import etomica.api.IBox;


/**
 * Event that conveys that an Atom has been added to a Box.
 */
public class BoxAtomAddedEvent extends BoxAtomEvent {

    public BoxAtomAddedEvent(IBox box, IAtom atom) {
        super(box, atom);
    }

    private static final long serialVersionUID = 1L;
}
