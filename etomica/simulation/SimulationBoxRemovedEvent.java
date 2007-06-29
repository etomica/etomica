package etomica.simulation;

import etomica.box.Box;

public class SimulationBoxRemovedEvent extends SimulationBoxEvent {

    public SimulationBoxRemovedEvent(Box box) {
        super(box);
    }

    private static final long serialVersionUID = 1L;
}
