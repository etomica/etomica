package etomica.simulation;

import etomica.api.IAtomType;

public class SimulationAtomTypeEvent extends SimulationEvent {

    public SimulationAtomTypeEvent(IAtomType atomType) {
        super();
        this.atomType = atomType;
    }

    public IAtomType getAtomType() {
        return atomType;
    }
    
    private final IAtomType atomType;
    private static final long serialVersionUID = 1L;
}
