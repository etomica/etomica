package etomica.box;


/**
 * Event that conveys that the maximum leaf index in a Box has changed
 * (or is about to change).
 */
public class BoxGlobalAtomLeafIndexEvent extends BoxEvent {

    public BoxGlobalAtomLeafIndexEvent(Box box, int maxIndex) {
        super(box);
        this.maxIndex = maxIndex;
    }
    
    public int getMaxIndex() {
        return maxIndex;
    }
    
    private final int maxIndex;
    private static final long serialVersionUID = 1L;
}
