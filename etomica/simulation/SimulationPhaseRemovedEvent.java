package etomica.simulation;

import etomica.phase.Phase;

public class SimulationPhaseRemovedEvent extends SimulationPhaseEvent {

    public SimulationPhaseRemovedEvent(Phase phase) {
        super(phase);
    }

    private static final long serialVersionUID = 1L;
}
