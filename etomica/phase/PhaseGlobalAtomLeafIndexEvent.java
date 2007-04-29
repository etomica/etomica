package etomica.phase;


/**
 * Event that conveys that the maximum leaf index in a Phase has changed
 * (or is about to change).
 */
public class PhaseGlobalAtomLeafIndexEvent extends PhaseEvent {

    public PhaseGlobalAtomLeafIndexEvent(Phase phase, int maxIndex) {
        super(phase);
        this.maxIndex = maxIndex;
    }
    
    public int getMaxIndex() {
        return maxIndex;
    }
    
    private final int maxIndex;
    private static final long serialVersionUID = 1L;
}
