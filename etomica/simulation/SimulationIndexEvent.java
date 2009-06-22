package etomica.simulation;

import etomica.api.ISimulation;
import etomica.api.ISimulationIndexEvent;

public class SimulationIndexEvent extends SimulationEvent implements ISimulationIndexEvent{

    private int index;
    
    public SimulationIndexEvent(ISimulation sim, int index) {
        super(sim);
        this.index = index;
    }
    
    public int getIndex() {
        return index;
    }
}
