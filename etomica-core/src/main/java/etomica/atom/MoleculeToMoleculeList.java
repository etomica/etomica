/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.atom;


/**
 * Interface for class that determines an AtomArrayList given an atom.
 */
public interface MoleculeToMoleculeList {

    /**
     * Returns the ArrayList that this instance associates with the given atom.
     * Should return null if the atom is null.
     */
    public IMoleculeList getMoleculeList(IMolecule atom);
}
