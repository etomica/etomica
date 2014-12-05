/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.api;

/**
 * Interface for a set of IAtoms.  The IAtomSet might contain 0, 1, 2 or many
 * IAtoms.
 * 
 * @author Andrew Schultz
 */
public interface IAtomList {

    /**
     * Returns the i-th atom, with numbering beginning from 0. 
     * If i is greater than count-1, throws an IllegalArgumentException.
     */
    public IAtom getAtom(int i);

    /**
     * @return the number of atoms in the set
     */
    public int getAtomCount();
}