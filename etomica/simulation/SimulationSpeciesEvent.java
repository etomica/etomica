package etomica.simulation;

import etomica.api.ISpecies;

public class SimulationSpeciesEvent extends SimulationEvent {

    public SimulationSpeciesEvent(ISpecies species) {
        super();
        this.species = species;
    }

    public ISpecies getSpecies() {
        return species;
    }
    
    private final ISpecies species;
    private static final long serialVersionUID = 1L;
}
