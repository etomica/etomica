/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.api;

/**
 * IAtomType identifies a set of atoms and defines properties of those atoms.
 * Properties include indicies used for tracking, mass and element.
 */
public interface IAtomType {

    /**
     * Informs the atom type what species contains the atom types.  This should
     * only be called by the species.
     */
    public void setSpecies(ISpecies newParent);

    /**
     * Informs the IAtomType what its index should be.  This should only be
     * called by the species.
     */
    public void setIndex(int newIndex);

    /**
     * Returns the index for this IAtomType, within the context of an
     * ISimulation.  The index is the IAtomType's position in the list of
     * atom types in the simulation.
     */
    public int getIndex();

    /**
     * Informs the atom type what its child index is.  This should only be called
     * by the species.
     */
    public void setChildIndex(int newChildIndex);

    /**
     * Returns the child index.  This is the index of the atom type within the
     * species.
     */
    public int getChildIndex();

    /**
     * Returns the species that contains the atom type
     */
    public ISpecies getSpecies();

    /**
     * Returns the value of the mass.
     */
    public double getMass();

    /**
     * Returns the reciprocal of the mass, 1.0/mass
     */
    public double rm();

    /**
     * Return the element for this atom type
     */
    public IElement getElement();

}