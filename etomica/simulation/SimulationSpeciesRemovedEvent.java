package etomica.simulation;

import etomica.api.ISpecies;

public class SimulationSpeciesRemovedEvent extends SimulationSpeciesEvent {

    public SimulationSpeciesRemovedEvent(ISpecies species) {
        super(species);
    }
    
    private static final long serialVersionUID = 1L;
}
