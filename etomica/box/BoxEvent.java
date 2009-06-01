package etomica.box;

import etomica.api.IBox;
import etomica.api.IBoxEvent;

/**
 * Event that conveys some happening with respect to a box or the things it contains.
 *
 * @see BoxListenerAdapter
 */
public class BoxEvent implements java.io.Serializable, IBoxEvent {
    
    public BoxEvent(IBox box) {
        this.box = box;
    }
    
    public IBox getBox() {
        return box;
    }
    
    protected IBox box;

    private static final long serialVersionUID = 1L;
}
    