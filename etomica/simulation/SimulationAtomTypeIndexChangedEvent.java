package etomica.simulation;

import etomica.api.IAtomType;

/**
 * Simulation event that notifies listeners that the index of an AtomType has
 * changed.
 * @author Andrew Schultz
 */
public class SimulationAtomTypeIndexChangedEvent extends SimulationAtomTypeEvent {

    SimulationAtomTypeIndexChangedEvent(IAtomType atomType, int oldIndex) {
        super(atomType);
        this.oldIndex = oldIndex;
    }
    
    /**
     * Returns the first AtomType removed
     */
    public int getOldIndex() {
        return oldIndex;
    }
    
    private final int oldIndex;
    private static final long serialVersionUID = 1L;
}
