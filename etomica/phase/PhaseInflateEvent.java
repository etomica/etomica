package etomica.phase;


/**
 * Event that conveys that the boundary dimensions of a phase have changed.
 */
public class PhaseInflateEvent extends PhaseEvent {

    public PhaseInflateEvent(Phase phase) {
        super(phase);
    }

    private static final long serialVersionUID = 1L;
}
