package etomica.simulation;

import etomica.species.Species;

public class SimulationSpeciesRemovedEvent extends SimulationSpeciesEvent {

    public SimulationSpeciesRemovedEvent(Species species) {
        super(species);
    }
    
    private static final long serialVersionUID = 1L;
}
