package etomica.simulation;

/**
 * Simulation event that notifies listeners that the number of Species in 
 * the simulation (and therefore the maximum index) has changed.
 * 
 * @author Rob Rassler
 */
public class SimulationSpeciesMaxIndexEvent extends SimulationEvent {

    SimulationSpeciesMaxIndexEvent(int maxIndex) {
        super();
        this.maxIndex = maxIndex;
    }
    
    /**
     * Returns the new maximum index
     */
    public int getMaxIndex() {
        return maxIndex;
    }
    
    private final int maxIndex;
    private static final long serialVersionUID = 1L;
}
