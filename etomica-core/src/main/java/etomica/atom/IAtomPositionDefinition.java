/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.atom;

import etomica.api.IMolecule;
import etomica.api.IVector;



/**
 * Returns a vector given an atom, thereby defining the position
 * of the atom or atom group.  Example implementations of this interface
 * are based on the center of mass, or on the position of the first
 * leaf atom in the group.
 */
public interface IAtomPositionDefinition {

    /**
     * Returns the defined position for the given atom, which 
     * may be an atom group.
     */
    public IVector position(IMolecule atom);
}
