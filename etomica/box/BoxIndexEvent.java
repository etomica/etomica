package etomica.box;

import etomica.api.IBox;
import etomica.api.IBoxIndexEvent;


/**
 * Event that conveys that the maximum leaf index in a Box has changed
 * (or is about to change).
 */
public class BoxIndexEvent extends BoxEvent implements IBoxIndexEvent {

    public BoxIndexEvent(IBox box, int _index) {
        super(box);
        this.index = _index;
    }

    public int getIndex() {
        return index;
    }
    
    protected int index = -1;
    private static final long serialVersionUID = 1L;
}
