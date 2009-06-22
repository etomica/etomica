package etomica.simulation;

import etomica.api.IAtomType;
import etomica.api.ISimulation;
import etomica.api.ISimulationAtomTypeEvent;

public class SimulationAtomTypeEvent extends SimulationEvent
                            implements ISimulationAtomTypeEvent {
    
    private final IAtomType atomType;
    private static final long serialVersionUID = 1L;
    
    public SimulationAtomTypeEvent(ISimulation sim, IAtomType atomType) {
        super(sim);
        this.atomType = atomType;
    }

    public IAtomType getAtomType() {
        return atomType;
    }

}
