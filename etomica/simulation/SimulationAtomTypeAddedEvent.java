package etomica.simulation;

import etomica.atom.AtomType;

public class SimulationAtomTypeAddedEvent extends SimulationAtomTypeEvent {

    public SimulationAtomTypeAddedEvent(AtomType newAtomType) {
        super(newAtomType);
    }

    private static final long serialVersionUID = 1L;
}
