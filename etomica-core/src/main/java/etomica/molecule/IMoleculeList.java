/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.molecule;

/**
 * Interface for a set of IMolecules.  The IMoleculeList might contain 0, 1, 2 or many
 * IMolecules.
 * 
 * @author Andrew Schultz
 */
public interface IMoleculeList {

    /**
     * Returns the i-th molecule, with numbering beginning from 0.
     * If i is greater than count-1, throws an IllegalArgumentException.
     */
    IMolecule get(int i);

    /**
     * @return the number of molecules in the list
     */
    int size();
}
