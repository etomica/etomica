package etomica.phase;

/**
 * Event that conveys some happening with respect to a phase or the things it contains.
 *
 * @see PhaseListener
 */
public class PhaseEvent implements java.io.Serializable {
    
    public PhaseEvent(Phase phase) {
        super();
        this.phase = phase;
    }
    
    public Phase getPhase() {
        return phase;
    }
    
    private final Phase phase;
    private static final long serialVersionUID = 1L;
}
    