package etomica.simulation;

import etomica.api.IAtomType;

public class SimulationAtomTypeAddedEvent extends SimulationAtomTypeEvent {

    public SimulationAtomTypeAddedEvent(IAtomType newAtomType) {
        super(newAtomType);
    }

    private static final long serialVersionUID = 1L;
}
