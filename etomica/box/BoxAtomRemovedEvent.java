package etomica.box;

import etomica.api.IAtom;
import etomica.api.IBox;
import etomica.api.IBoxAtomRemovedEvent;


/**
 * Event that conveys that an Atom has been removed from a Box.
 */
public class BoxAtomRemovedEvent extends BoxAtomEvent implements IBoxAtomRemovedEvent {

    public BoxAtomRemovedEvent(IBox box, IAtom atom) {
        super(box, atom);
    }

    private static final long serialVersionUID = 1L;
}
