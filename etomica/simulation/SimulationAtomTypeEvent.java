package etomica.simulation;

import etomica.api.IAtomTypeLeaf;

public class SimulationAtomTypeEvent extends SimulationEvent {

    public SimulationAtomTypeEvent(IAtomTypeLeaf atomType) {
        super();
        this.atomType = atomType;
    }

    public IAtomTypeLeaf getAtomType() {
        return atomType;
    }
    
    private final IAtomTypeLeaf atomType;
    private static final long serialVersionUID = 1L;
}
