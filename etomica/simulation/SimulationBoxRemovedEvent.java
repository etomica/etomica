package etomica.simulation;

import etomica.api.IBox;

public class SimulationBoxRemovedEvent extends SimulationBoxEvent {

    public SimulationBoxRemovedEvent(IBox box) {
        super(box);
    }

    private static final long serialVersionUID = 1L;
}
