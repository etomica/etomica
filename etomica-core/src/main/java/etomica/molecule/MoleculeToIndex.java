/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.molecule;


/**
 * Interface for class that associates an integer with an atom.
 */
public interface MoleculeToIndex {

    /**
     * Returns an integer that this instance associates with the given atom.
     * Should return -1 if an appropriate index can not be determined.
     */
    public int getIndex(IMolecule atom);
}
