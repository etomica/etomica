package etomica.atom;

import etomica.simulation.SimulationEvent;

class SimulationAtomTypeCompactedEvent extends SimulationEvent {

    SimulationAtomTypeCompactedEvent(int startIndex, int stopIndex) {
        this.startIndex = startIndex;
        this.stopIndex = stopIndex;
    }
    
    /**
     * Returns the first AtomType removed
     */
    int getStartIndex() {
        return startIndex;
    }
    
    /**
     * Returns the first AtomType not removed
     */
    int getStopIndex() {
        return stopIndex;
    }

    private final int startIndex, stopIndex;
    private static final long serialVersionUID = 1L;
}
