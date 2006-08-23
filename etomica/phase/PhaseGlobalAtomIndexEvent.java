package etomica.phase;


/**
 * Event that conveys that the maximum global index in a Phase has changed.
 */
public class PhaseGlobalAtomIndexEvent extends PhaseEvent {

    public PhaseGlobalAtomIndexEvent(Phase phase, int maxIndex) {
        super(phase);
        this.maxIndex = maxIndex;
    }
    
    public int getMaxIndex() {
        return maxIndex;
    }
    
    private final int maxIndex;
    private static final long serialVersionUID = 1L;
}
