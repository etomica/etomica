/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.molecule;

import etomica.box.Box;
import etomica.space.IOrientation;
import etomica.space.Space;
import etomica.space.Vector;

/**
 * Interface for a class that can calculate the orientation of a molecule or
 * to set the orientation of a molecule.
 *
 * @author Andrew Schultz
 */
public interface OrientationCalc {

    IOrientation makeOrientation(Space space);

    int getDOF(Space space);

    Vector getMomentOfInertia(IMolecule molecule);

    /**
     * Calculates the orientation of the given molecule and stores that
     * in the given orientation.
     */
    void calcOrientation(IMolecule molecule, IOrientation orientation);

    /**
     * Sets the orientation of the given molecule to be equal to the given
     * orientation.  This typically involves changing the position of the atoms
     * without changing the molecules position.
     */
    void setOrientation(IMolecule molecule, Box box, IOrientation orientation);

    Vector getAngularMomentum(IMolecule molecule, Vector com, Box box);

    void setAngularMomentum(IMolecule molecule, Vector com, Box box, Vector L);

    Vector angularMomentumToVelocity(IMolecule molecule, Space space, Vector L);

    Vector angularMomentumToVelocity(IOrientation orientation, IMolecule molecule, Space space, Vector L);

    Vector angularVelocityToMomentum(IMolecule molecule, Box box, Vector omega);

    Vector bodyToSpace(IMolecule molecule, Space space, Vector v);

    Vector spaceToBody(IMolecule molecule, Space space, Vector v);
}