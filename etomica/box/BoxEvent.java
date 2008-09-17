package etomica.box;

import etomica.api.IBox;
import etomica.api.IEvent;

/**
 * Event that conveys some happening with respect to a box or the things it contains.
 *
 * @see BoxListener
 */
public class BoxEvent implements IEvent, java.io.Serializable {
    
    public BoxEvent(IBox box) {
        super();
        this.box = box;
    }
    
    public IBox getBox() {
        return box;
    }
    
    private final IBox box;
    private static final long serialVersionUID = 1L;
}
    