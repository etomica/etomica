package etomica.simulation;

import etomica.species.Species;

public class SimulationSpeciesAddedEvent extends SimulationSpeciesEvent {

    public SimulationSpeciesAddedEvent(Species species) {
        super(species);
    }
    
    private static final long serialVersionUID = 1L;
}
