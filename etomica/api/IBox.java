package etomica.api;

public interface IBox {

    /**
     * Informs the IBox what its index is.  This should only be called by the
     * ISimulation.
     */
	public void setIndex(int newIndex);

	/**
	 * Returns the IBox's index.
	 */
	public int getIndex();

	/**
	 * Adds a new molecule of the 
	 * @param species
	 * @return
	 */
	public IMolecule addNewMolecule(ISpecies species);

	public void addMolecule(IMolecule molecule);

	public void removeMolecule(IMolecule molecule);

	public void setNMolecules(ISpecies species, int n);

	public int getNMolecules(ISpecies species);

	public IAtomSet getMoleculeList(ISpecies species);

	public IAtomSet getMoleculeList();

	public void setBoundary(IBoundary b);

	public IBoundary getBoundary();

	public void setDimensions(IVector d);

	public void setDensity(double rho);

	public IBoxEventManager getEventManager();

	public void addSpeciesNotify(ISpecies species);

	/**
	 * Notifies the SpeciesMaster that a Species has been removed.  This method
	 * should only be called by the SpeciesManager.
	 */
	public void removeSpeciesNotify(ISpecies species);

	public IAtomSet getLeafList();

	public int requestGlobalIndex();

	public int getMaxGlobalIndex();

	public void addAtomNotify(IAtom newAtom);

	//updating of leaf atomList may not be efficient enough for repeated
	// use, but is probably ok
	public void removeAtomNotify(IAtom oldAtom);

	public int getLeafIndex(IAtom atomLeaf);

}