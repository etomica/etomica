package etomica.api;

public interface IBox {


	public abstract void resetIndex(ISimulation sim);

	public abstract int getIndex();

	public abstract IMolecule addNewMolecule(ISpecies species);

	public abstract void addMolecule(IMolecule molecule);

	public abstract void removeMolecule(IMolecule molecule);

	public abstract void setNMolecules(ISpecies species, int n);

	public abstract int getNMolecules(ISpecies species);

	public abstract IAtomSet getMoleculeList(ISpecies species);

	public abstract IAtomSet getMoleculeList();

	public abstract void setBoundary(IBoundary b);

	public abstract IBoundary getBoundary();

	public abstract void setDimensions(IVector d);

	public abstract void setDensity(double rho);

	public abstract IBoxEventManager getEventManager();

	public abstract void addSpeciesNotify(ISpecies species);

	/**
	 * Notifies the SpeciesMaster that a Species has been removed.  This method
	 * should only be called by the SpeciesManager.
	 */
	public abstract void removeSpeciesNotify(ISpecies species);

	public abstract IAtomSet getLeafList();

	public abstract int requestGlobalIndex();

	public abstract int getMaxGlobalIndex();

	public abstract void addAtomNotify(IAtom newAtom);

	//updating of leaf atomList may not be efficient enough for repeated
	// use, but is probably ok
	public abstract void removeAtomNotify(IAtom oldAtom);

	public abstract int getLeafIndex(IAtom atomLeaf);

}