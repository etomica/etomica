package etomica.api;


public interface ISpeciesManager {

    /**
     * Adds species to the list of all species in the simulation, and
     * adds new species agent to every box currently in simulation.
     * This is called by the Species constructor.
     * 
     * @return the index assigned to the new species
     */
    public void addSpecies(ISpecies species);

    public void boxAddedNotify(IBox newBox);

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