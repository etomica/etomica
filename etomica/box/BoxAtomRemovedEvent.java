package etomica.box;

import etomica.atom.IAtom;


/**
 * Event that conveys that an Atom has been removed from a Box.
 */
public class BoxAtomRemovedEvent extends BoxAtomEvent {

    public BoxAtomRemovedEvent(Box box, IAtom atom) {
        super(box, atom);
    }

    private static final long serialVersionUID = 1L;
}
