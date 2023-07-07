/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.molecule;


import etomica.atom.IAtomList;
import etomica.species.ISpecies;
import etomica.util.Statefull;

import java.util.ArrayList;

/**
 * Interface for an IMolecule
 */
public interface IMolecule extends Statefull {

    /**
     * Returns this IMolecule's index, which is its place in the list of
     * molecules of its ISpecies in its Box.
     */
    int getIndex();

    /**
     * Informs the IMolecule of its index.
     */
    void setIndex(int index);

    /**
     * Returns the atoms in the molecule as an IAtomList.
     */
    IAtomList getChildList();

    /**
     * Returns the ISpecies of this IMolecule.
     */
    ISpecies getType();

    void copyCoordinatesFrom(IMolecule molecule);

}
