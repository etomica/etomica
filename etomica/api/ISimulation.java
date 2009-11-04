package etomica.api;

public interface ISimulation {

    /**
     * Adds a Box to the simulation.  This method should not be called if
     * newBox is already held by the simulation.
     */
    public void addBox(IBox newBox);

    /**
     * Removes a Box to the simulation.  This method should not be called if
     * oldBox is not held by the simulation.
     */
    public void removeBox(IBox oldBox);

    /**
     * Returns number of boxes contained in the Simulation
     */
    public int getBoxCount();

    /**
     * Returns Box specified by index contained in the Simulation
     */
    public IBox getBox(int index);

    /**
     * Returns the Simulation's random number generator.
     */
    public IRandom getRandom();

    /**
     * Returns the Simulation's event manager, which fires events for
     * Boxs and Species being added and removed.
     */
    public ISimulationEventManager getEventManager();

    /**
     * Adds species to the list of all species in the simulation, and
     * adds new species agent to every box currently in simulation.
     * This is called by the Species constructor.
     * 
     * @return the index assigned to the new species
     */
    public void addSpecies(ISpecies species);

    /**
     * Removes the given AtomTypes associated with the given Species from the 
     * Simulation and does cleanup, including renumbering indices and firing 
     * AtomType-related event notifications.
     */
    public void removeSpecies(ISpecies removedSpecies);

    /**
     * Returns the number of Species in the Simulation.
     */
    public int getSpeciesCount();

    /**
     * Returns the Species in the Simulation for the specified index.
     */
    public ISpecies getSpecies(int index);
}
