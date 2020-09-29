/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.molecule;


import etomica.atom.IAtomList;
import etomica.species.ISpecies;

/**
 * Interface for an IMolecule
 */
public interface IMolecule {

    /**
     * Returns this IMolecule's index, which is its place in the list of
     * molecules of its ISpecies in its Box.
     */
    public int getIndex();

    /**
     * Informs the IMolecule of its index.
     */
    public void setIndex(int index);

    void setGlobalIndex(int idx);

    int getGlobalIndex();

    /**
     * Returns the atoms in the molecule as an IAtomList.
     */
    public IAtomList getChildList();

    /**
     * Returns the ISpecies of this IMolecule.
     */
    public ISpecies getType();

    void copyFrom(IMolecule other);
}
