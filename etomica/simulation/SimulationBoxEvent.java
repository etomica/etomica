package etomica.simulation;

import etomica.api.IBox;
import etomica.api.ISimulation;
import etomica.api.ISimulationBoxEvent;

public class SimulationBoxEvent extends SimulationEvent implements ISimulationBoxEvent {

    private final IBox box;
    private static final long serialVersionUID = 1L;
    
    public SimulationBoxEvent(ISimulation sim, IBox box) {
        super(sim);
        this.box = box;
    }

    public IBox getBox() {
        return box;
    }

}
