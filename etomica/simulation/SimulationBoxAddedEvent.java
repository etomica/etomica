package etomica.simulation;

import etomica.api.IBox;

public class SimulationBoxAddedEvent extends SimulationBoxEvent {

    public SimulationBoxAddedEvent(IBox box) {
        super(box);
    }

    private static final long serialVersionUID = 1L;
}
