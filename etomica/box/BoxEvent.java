package etomica.box;

import etomica.api.IBox;
import etomica.api.IBoxEvent;
import etomica.api.IEvent;

/**
 * Event that conveys some happening with respect to a box or the things it contains.
 *
 * @see BoxListener
 */
public class BoxEvent implements IEvent, java.io.Serializable, IBoxEvent {
    
    public BoxEvent(IBox box) {
        super();
        this.box = box;
    }
    
    /* (non-Javadoc)
     * @see etomica.box.IBoxEvent#getBox()
     */
    public IBox getBox() {
        return box;
    }
    
    private final IBox box;
    private static final long serialVersionUID = 1L;
}
    