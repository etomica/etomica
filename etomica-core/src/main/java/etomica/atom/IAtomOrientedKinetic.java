/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.atom;

import etomica.space.Vector;

/**
 * Interface for an Atom that has a position, orientation, velocity and angular
 * velocity.
 */
public interface IAtomOrientedKinetic extends IAtomKinetic, IAtomOriented {

    //XXX angular velocity is not a vector.  enjoy!
    public Vector getAngularVelocity(); //angular velocity vector in space-fixed frame
    
}
