/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.atom;

import etomica.space.Vector;



/**
 * Returns a vector given a molecule, thereby defining the position
 * of the molecule.  Example implementations of this interface
 * are based on the center of mass, or on the position of the first
 * child atom in the molecule.
 */
public interface IMoleculePositionDefinition {

    /**
     * Returns the defined position for the given molecule.
     * @param molecule a group of atoms for which to define a position
     */
    Vector position(IMolecule molecule);
}
