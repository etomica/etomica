/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.atom;

import java.util.List;

/**
 * Interface for a list of IAtoms.  The IAtomList might contain 0, 1, 2 or many
 * IAtoms.
 * 
 * @author Andrew Schultz
 */
public interface IAtomList extends List<IAtom> {

    /**
     * Returns the i-th atom, with numbering beginning from 0.
     * @param i the index of the desired atom in the list
     *          if i is greater than count-1, throws an IllegalArgumentException.
     * @return the requested atom
     */
    IAtom get(int i);

    /**
     * @return the number of atoms in the list
     */
    int size();

}
