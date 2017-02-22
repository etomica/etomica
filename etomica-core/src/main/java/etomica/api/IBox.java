/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.api;

public interface IBox {

    /**
     * Informs the IBox what its index is.  This should only be called by the
     * ISimulation.
     */
    public void setIndex(int newIndex);

    /**
     * Returns the IBox's index.  The index corresponds to the box's position
     * in the simulation's list of IBoxes.  The index of the first IBox is 0.
     * The index of the last IBox is n-1, where n is the number of IBoxes.
     */
    public int getIndex();

    /**
     * Adds the given molecule to the this box.  The molecule should not
     * already be in this box and should not be in another IBox.  The molecule
     * should be a member of an ISpecies which has been added to the
     * ISimulation.
     */
    public void addMolecule(IMolecule molecule);

    /**
     * Removes the given molecule from this box.  The molecule must be held
     * by the box before this method is called.
     */
    public void removeMolecule(IMolecule molecule);

    /**
     * Sets the number of molecules in this box of the given ISpecies to n.
     * Molecules are added to or removed from the box to achieve the desired
     * number.
     */
    public void setNMolecules(ISpecies species, int n);

    /**
     * Returns the number of molecules in this box of the given ISpecies.
     */
    public int getNMolecules(ISpecies species);

    /**
     * Returns the list of molecules of the given species as an IAtomSet that
     * are in this box.
     */
    public IMoleculeList getMoleculeList(ISpecies species);

    /**
     * Returns a list of all molecules in this box as an IAtomSet.
     */
    public IMoleculeList getMoleculeList();

    /**
     * Returns the list of atoms contained in this box.
     */
    public IAtomList getLeafList();

    /**
     * Sets the box's boundary to the given IBoundary.
     */
    public void setBoundary(IBoundary newBoundary);

    /**
     * Returns the box's boundary.
     */
    public IBoundary getBoundary();

    /**
     * Returns the event manager for this box.  
     * @return
     */
    public IBoxEventManager getEventManager();

    /**
     * Notifies the IBox that the given species has been added to the
     * simulation.  This method should only be called by the simulation.
     */
    public void addSpeciesNotify(ISpecies species);

    /**
     * Notifies the IBox that a Species has been removed.  This method should
     * only be called by the simulation.  This triggers the removal of all
     * molecules of the given species from this box.
     */
    public void removeSpeciesNotify(ISpecies species);
}
