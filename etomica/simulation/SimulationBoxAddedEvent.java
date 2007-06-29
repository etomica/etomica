package etomica.simulation;

import etomica.box.Box;

public class SimulationBoxAddedEvent extends SimulationBoxEvent {

    public SimulationBoxAddedEvent(Box box) {
        super(box);
    }

    private static final long serialVersionUID = 1L;
}
