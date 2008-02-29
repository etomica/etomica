package etomica.box;

import etomica.api.IBox;


/**
 * Event that conveys that the boundary dimensions of a box have changed.
 */
public class BoxInflateEvent extends BoxEvent {

    public BoxInflateEvent(IBox box) {
        super(box);
    }

    private static final long serialVersionUID = 1L;
}
