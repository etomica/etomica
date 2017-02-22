/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.api;


public interface IAtomType {

    public void setSpecies(ISpecies newParent);

    /**
     * Informs the ISpecies what its index should be.  This should only be
     * called by the SpeciesManager.
     */
    public void setIndex(int newIndex);

    /**
     * Returns the index for this IAtomType, within the context of an
     * ISimulation.  The index is the IAtomType's position in the list of
     * atom types in the simulation.
     */
    public int getIndex();
    
    public void setChildIndex(int newChildIndex);

    public int getChildIndex();

    public ISpecies getSpecies();

    /**
     * Returns the value of the mass.
     */
    public double getMass();

    /**
     * Returns the reciprocal of the mass, 1.0/mass
     */
    public double rm();

    public IElement getElement();

}