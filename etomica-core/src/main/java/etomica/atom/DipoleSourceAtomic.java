/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.atom;

import etomica.space.Vector;

/**
 * Interface for something that can calculate the dipole of an atom.
 */
public interface DipoleSourceAtomic {
    
    /**
     * Returns the dipole of the given atom.
     * The method is likely to through an exception if a atom of the wrong
     * type is passed in.
     */
    Vector getDipole(IAtom atom);
}
