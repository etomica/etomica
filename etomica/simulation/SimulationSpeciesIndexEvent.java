package etomica.simulation;

import etomica.api.ISimulation;
import etomica.api.ISimulationSpeciesIndexEvent;
import etomica.api.ISpecies;

public class SimulationSpeciesIndexEvent extends SimulationSpeciesEvent
                                 implements ISimulationSpeciesIndexEvent {

    private int index;
    
    public SimulationSpeciesIndexEvent(ISimulation sim, ISpecies species,
                                       int index) {
        super(sim, species);
        this.index = index;
    }
    
    public int getIndex() {
        return index;
    }
}
