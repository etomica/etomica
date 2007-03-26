package etomica.simulation;


/**
 * Simulation event that notifies listeners that the number of AtomTypes in 
 * the simulation (and therefore the maximum index) has changed.
 * 
 * @author Andrew Schultz
 */
public class SimulationAtomTypeMaxIndexEvent extends SimulationEvent {

    SimulationAtomTypeMaxIndexEvent(int maxIndex) {
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
