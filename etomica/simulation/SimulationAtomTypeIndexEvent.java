package etomica.simulation;

import etomica.api.IAtomType;
import etomica.api.ISimulation;
import etomica.api.ISimulationAtomTypeIndexEvent;

public class SimulationAtomTypeIndexEvent extends SimulationAtomTypeEvent
                           implements ISimulationAtomTypeIndexEvent {

    private int index;
    
    public SimulationAtomTypeIndexEvent(ISimulation sim, IAtomType at, int index) {
        super(sim, at);
        this.index = index;
    }
    
    public int getIndex() {
        return index;
    }
}
