package etomica.simulation;

import etomica.api.ISimulation;
import etomica.api.ISimulationSpeciesEvent;
import etomica.api.ISpecies;

public class SimulationSpeciesEvent extends SimulationEvent
                                   implements ISimulationSpeciesEvent {

    private final ISpecies species;
    private static final long serialVersionUID = 1L;
    
    public SimulationSpeciesEvent(ISimulation sim, ISpecies species) {
        super(sim);
        this.species = species;
    }

    public ISpecies getSpecies() {
        return species;
    }

}
