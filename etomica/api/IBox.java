package etomica.api;

import etomica.atom.AtomSet;
import etomica.atom.IAtom;
import etomica.atom.IMolecule;
import etomica.box.BoxEventManager;
import etomica.simulation.ISimulation;
import etomica.space.Boundary;
import etomica.species.ISpecies;

public interface IBox {

	/**
	 * Resets the Box's index.  This should only need to be called from the
	 * Simulation class.
	 * 
	 * @param sim  The Simulation to which this Box was added.  Passing null
	 *             notifies the Box that it was removed from the Simulation
	 *             (the index is set to 0).
	 */
	public abstract void resetIndex(ISimulation sim);

	public abstract int getIndex();

	public abstract IMolecule addNewMolecule(ISpecies species);

	public abstract void addMolecule(IMolecule molecule);

	public abstract void removeMolecule(IMolecule molecule);

	/**
	 * Sets the number of molecules for this species.  Molecules are either
	 * added or removed until the given number is obtained.  Takes no action
	 * at all if the new number of molecules equals the existing number.
	 *
	 * @param n  the new number of molecules for this species
	 */
	public abstract void setNMolecules(ISpecies species, int n);

	public abstract int getNMolecules(ISpecies species);

	public abstract AtomSet getMoleculeList(ISpecies species);

	public abstract AtomSet getMoleculeList();

	/**
	 * Sets the boundary object of the box.
	 */
	public abstract void setBoundary(Boundary b);

	/**
	 * Returns the current boundary instance.
	 * 
	 * @return The current instance of the boundary class
	 */
	public abstract Boundary getBoundary();

	public abstract BoxEventManager getEventManager();

	/**
	 * Returns an AtomArrayList containing the leaf atoms in the Box
	 */
	public abstract AtomSet getLeafList();

	/**
	 * Returns a "global" index for the Box.  This method should only be
	 * called by Atom.
	 */
	public abstract int requestGlobalIndex();

	/**
	 * Returns the maximum global index for the Box.
	 */
	public abstract int getMaxGlobalIndex();

	/**
	 * Returns the index of the given leaf atom within the SpeciesMaster's
	 * leaf list.  The given leaf atom must be in the SpeciesMaster's Box. 
	 */
	public abstract int getLeafIndex(IAtom atomLeaf);

}