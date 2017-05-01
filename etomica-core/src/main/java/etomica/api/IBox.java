/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.api;

/**
 * IBox holds a boundary and a list of molecules.  The boundary defines the
 * spatial extent of the box and conditions at the edges of the boundary.  The
 * molecules in the box interact via potentials, but do not interact with
 * molecules in other boxes.  IBox contains methods to add and remove molecules
 * and an event manager that notifies listeners about those events.
 */
public interface IBox {

    /**
     * Informs the IBox what its index is.  This should only be called by the
     * ISimulation.
     *
     * @param newIndex the box's new index
     */
    public void setIndex(int newIndex);

    /**
     * @return the IBox's index.  The index corresponds to the box's position
     * in the simulation's list of IBoxes.  The index of the first IBox is 0.
     * The index of the last IBox is n-1, where n is the number of IBoxes.
     */
    public int getIndex();

    /**
     * Adds the given molecule to the this box.  The molecule should not
     * already be in this box and should not be in another IBox.  The molecule
     * should be a member of an ISpecies which has been added to the
     * ISimulation.
     *
     * @param molecule the molecule to be added to the box
     */
    public void addMolecule(IMolecule molecule);

    /**
     * Removes the given molecule from this box.  The molecule must be held
     * by the box before this method is called.
     *
     * @param molecule the molecule to be removed from the box
     */
    public void removeMolecule(IMolecule molecule);

    /**
     * Sets the number of molecules in this box of the given ISpecies to n.
     * Molecules are added to or removed from the box to achieve the desired
     * number.
     *
     * @param species the species whose number of molecules should be changed
     * @param n the desired number of molecules
     */
    public void setNMolecules(ISpecies species, int n);

    /**
     * @return the number of molecules in this box of the given ISpecies.
     */
    public int getNMolecules(ISpecies species);

    /**
     * Returns the list of molecules of the given species.
     *
     * @param species the species
     * @return the requested list of molecules
     */
    public IMoleculeList getMoleculeList(ISpecies species);

    /**
     * @return a list of all molecules in this box.
     */
    public IMoleculeList getMoleculeList();

    /**
     * @return the list of atoms contained in this box.
     */
    public IAtomList getLeafList();

    /**
     * Sets the box's boundary to the given IBoundary.
     *
     * @param newBoundary the new boundary
     */
    public void setBoundary(IBoundary newBoundary);

    /**
     * @return the box's boundary.
     */
    public IBoundary getBoundary();

    /**
     * @return the event manager for this box.
     */
    public IBoxEventManager getEventManager();

    /**
     * Notifies the IBox that the given species has been added to the
     * simulation.  This method should only be called by the simulation.
     *
     * @param species the added species
     */
    public void addSpeciesNotify(ISpecies species);

    /**
     * Notifies the IBox that a Species has been removed.  This method should
     * only be called by the simulation.  This triggers the removal of all
     * molecules of the given species from this box.
     *
     * @param species the removed species
     */
    public void removeSpeciesNotify(ISpecies species);
}
