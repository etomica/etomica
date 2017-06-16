/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.molecule;

import etomica.box.Box;

/**
 * Interface for objects when return atoms (meeting some specification)
 * from a box.
 */
public interface MoleculeSource {
    
    /**
     * sets the Box the source should pull Atoms from.
     * Box should not be null
     */
    public void setBox(Box p);

    /**
     * Returns an atom.  Will return null if there are no appropriate atoms in 
     * the given box.
     */
    public IMolecule getMolecule();
}
