package etomica.box;

import etomica.api.IAtom;
import etomica.api.IBox;
import etomica.api.IBoxAtomIndexEvent;


/**
 * Event that conveys that the maximum global index in a Box has changed.
 */
public class BoxAtomIndexEvent extends BoxAtomEvent implements IBoxAtomIndexEvent {

    public BoxAtomIndexEvent(IBox box, IAtom atom, int _index) {
        super(box, atom);
        this.index = _index;
    }

    public int getIndex() {
        return index;
    }
    
    protected int index = -1;
    private static final long serialVersionUID = 1L;
}
