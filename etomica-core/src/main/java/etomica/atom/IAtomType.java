/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.atom;

import etomica.api.IElement;
import etomica.api.ISpecies;

/**
 * IAtomType identifies a set of atoms and defines properties of those atoms.
 * Properties include indices used for tracking, mass and element.
 */
public interface IAtomType {

    /**
     * Informs the atom type what species contains the atom types.  This should
     * only be called by the species.
     * @param newParent the atom type's new parent
     */
    public void setSpecies(ISpecies newParent);

    /**
     * Informs the IAtomType what its index should be.  This should only be
     * called by the species.
     * @param newIndex the atom type's new index
     */
    public void setIndex(int newIndex);

    /**
     * @return the index for this IAtomType, within the context of an
     * ISimulation.  The index is the IAtomType's position in the list of
     * atom types in the simulation.
     */
    public int getIndex();

    /**
     * Informs the atom type what its child index is.  This should only be called
     * by the species.
     * @param newChildIndex the atom type's new child index
     */
    public void setChildIndex(int newChildIndex);

    /**
     * @return the child index.  This is the index of the atom type within the
     * species.
     */
    public int getChildIndex();

    /**
     * @return the species that contains the atom type
     */
    public ISpecies getSpecies();

    /**
     * @return the value of the mass.
     */
    public double getMass();

    /**
     * @return the reciprocal of the mass, 1.0/mass
     */
    public double rm();

    /**
     * @return the element for this atom type
     */
    public IElement getElement();

}
