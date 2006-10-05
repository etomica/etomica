package etomica.atom;

import etomica.simulation.SimulationAtomTypeEvent;

/**
 * Simulation event that notifies listeners that the index of an AtomType has
 * changed.
 * @author Andrew Schultz
 */
class SimulationAtomTypeIndexChangedEvent extends SimulationAtomTypeEvent {

    SimulationAtomTypeIndexChangedEvent(AtomType atomType, int oldIndex) {
        super(atomType);
        this.oldIndex = oldIndex;
    }
    
    /**
     * Returns the first AtomType removed
     */
    int getOldIndex() {
        return oldIndex;
    }
    
    private final int oldIndex;
    private static final long serialVersionUID = 1L;
}
