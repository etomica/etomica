package etomica.simulation;

import etomica.phase.Phase;

public class SimulationPhaseAddedEvent extends SimulationPhaseEvent {

    public SimulationPhaseAddedEvent(Phase phase) {
        super(phase);
    }

    private static final long serialVersionUID = 1L;
}
