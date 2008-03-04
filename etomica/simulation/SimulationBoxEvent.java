package etomica.simulation;

import etomica.api.IBox;

public class SimulationBoxEvent extends SimulationEvent {

    public SimulationBoxEvent(IBox box) {
        super();
        this.box = box;
    }

    public IBox getBox() {
        return box;
    }
    
    private final IBox box;
    private static final long serialVersionUID = 1L;
}
