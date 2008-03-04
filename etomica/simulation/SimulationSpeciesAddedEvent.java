package etomica.simulation;

import etomica.api.ISpecies;

public class SimulationSpeciesAddedEvent extends SimulationSpeciesEvent {

    public SimulationSpeciesAddedEvent(ISpecies species) {
        super(species);
    }
    
    private static final long serialVersionUID = 1L;
}
