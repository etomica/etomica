package etomica.simulation;

import etomica.box.Box;

public class SimulationBoxEvent extends SimulationEvent {

    public SimulationBoxEvent(Box box) {
        super();
        this.box = box;
    }

    public Box getBox() {
        return box;
    }
    
    private final Box box;
    private static final long serialVersionUID = 1L;
}
