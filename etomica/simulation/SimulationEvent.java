package etomica.simulation;

import etomica.api.ISimulation;
import etomica.api.ISimulationEvent;


public class SimulationEvent implements ISimulationEvent, java.io.Serializable {
    
    private ISimulation simulation;
    private static final long serialVersionUID = 1L;
    
    public SimulationEvent(ISimulation sim) {
    	simulation = sim;
    }
    
    public ISimulation getSimulation() {
        return simulation;
    }
}