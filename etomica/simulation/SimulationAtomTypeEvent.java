package etomica.simulation;

import etomica.atom.AtomType;

public class SimulationAtomTypeEvent extends SimulationEvent {

    public SimulationAtomTypeEvent(AtomType atomType) {
        super();
        this.atomType = atomType;
    }

    public AtomType getAtomType() {
        return atomType;
    }
    
    private final AtomType atomType;
    private static final long serialVersionUID = 1L;
}
