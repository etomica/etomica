package etomica.api;


public interface ISpeciesManager {

	/**
	 * Adds species to the list of all species in the simulation, and
	 * adds new species agent to every box currently in simulation.
	 * This is called by the Species constructor.
	 * 
	 * @return the index assigned to the new species
	 */
	public abstract void addSpecies(ISpecies species);

	public abstract void boxAddedNotify(IBox newBox);

	/**
	 * Removes the given AtomTypes associated with the given Species from the 
	 * Simulation and does cleanup, including renumbering indices and firing 
	 * AtomType-related event notifications.
	 */
	public abstract boolean removeSpecies(ISpecies removedSpecies);

	/**
	 * Returns an array of the Species in the Simulation.
	 */
	public abstract ISpecies[] getSpecies();

	/**
	 * This method notifies the SpeciesManager that the give atomType was added
	 * to the system.  This method should be called by the AtomType at the top
	 * of the AtomType hierarchy whenver it receives notification of a new
	 * AtomType.
	 */
	public abstract void atomTypeAddedNotify(IAtomType newChildType);

	public abstract int requestTypeIndex();

	public abstract void atomTypeRemovedNotify(IAtomType removedType);

	public abstract ISpecies[] getMoleculeTypes();

}