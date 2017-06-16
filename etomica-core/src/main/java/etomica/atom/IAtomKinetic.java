/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.atom;


import etomica.space.Vector;

/**
 * Interface for an atom that holds vectors for velocity.
 */
public interface IAtomKinetic extends IAtom {

    /**
     * @return the velocity of the IAtom.  Modifying the returned vector will
     * alter the IAtom's velocity.
     */
    public Vector getVelocity();
 
}
