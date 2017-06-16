/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.molecule;

import etomica.space.Vector;

/**
 * Interface for something that can calculate the dipole of a molecule.
 */
public interface DipoleSource {
    
    /**
     * Returns the dipole of the given molecule.
     * The method is likely to through an exception if a molecule of the wrong
     * type is passed in.
     */
    public Vector getDipole(IMolecule molecule);
}
