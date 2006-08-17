package etomica.simulation;

import etomica.phase.Phase;

public class SimulationPhaseEvent extends SimulationEvent {

    public SimulationPhaseEvent(Phase phase) {
        super();
        this.phase = phase;
    }

    public Phase getPhase() {
        return phase;
    }
    
    private final Phase phase;
    private static final long serialVersionUID = 1L;
}
