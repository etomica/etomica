/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.molecule;

import etomica.box.Box;
import etomica.space.Vector;
import etomica.space3d.IOrientationFull3D;

/**
 * Interface for a class that can calculate the orientation of a molecule or
 * to set the orientation of a molecule.
 *
 * @author Andrew Schultz
 */
public interface OrientationCalc {

    /**
     * Calculates the orientation of the given molecule and stores that
     * in the given orientation.
     */
    void calcOrientation(IMolecule molecule, IOrientationFull3D orientation);

    /**
     * Sets the orientation of the given molecule to be equal to the given
     * orientation.  This typically involves changing the position of the atoms
     * without changing the molecules position.
     */
    void setOrientation(IMolecule molecule, Box box, IOrientationFull3D orientation);

    Vector getAngularMomentum(IMolecule molecule, Vector com, Box box);

    void setAngularMomentum(IMolecule molecule, Vector com, Box box, Vector L);
}