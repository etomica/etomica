package etomica.simulation;

import etomica.api.ISpecies;

public class SimulationSpeciesIndexChangedEvent extends SimulationSpeciesEvent {

    SimulationSpeciesIndexChangedEvent(ISpecies species, int oldIndex) {
        super(species);
        this.oldIndex = oldIndex;
    }
    
    /**
     * Returns the first Species removed
     */
    public int getOldIndex() {
        return oldIndex;
    }
    
    private final int oldIndex;
    private static final long serialVersionUID = 1L;
}
