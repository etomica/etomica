package etomica.simulation;

import etomica.species.Species;

public class SimulationSpeciesEvent extends SimulationEvent {

    public SimulationSpeciesEvent(Species species) {
        super();
        this.species = species;
    }

    public Species getSpecies() {
        return species;
    }
    
    private final Species species;
    private static final long serialVersionUID = 1L;
}
