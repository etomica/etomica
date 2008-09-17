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
    public IEventManager getEventManager();

    /**
     * Returns the SpeciesManager, which tracks the Species in the Simulation.
     */
    public ISpeciesManager getSpeciesManager();

    /**
     * @return Returns a flag indicating whether the simulation involves molecular dynamics.
     */
    public boolean isDynamic();
}