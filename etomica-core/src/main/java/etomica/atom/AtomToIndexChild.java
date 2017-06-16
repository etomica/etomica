/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.atom;

import java.io.Serializable;


/**
 * Defines the index as the Atom's node's index.
 * @author andrew
 */
public class AtomToIndexChild implements AtomToIndex, Serializable {

    /**
     * @throws NullPointerException if the atom is null.
     */
    public int getIndex(IAtom atom) {
        return atom.getIndex();
    }
    
    private static final long serialVersionUID = 1L;

}
