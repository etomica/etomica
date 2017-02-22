/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.api;


/**
 * Interface for an atom that holds vectors for velocity.
 */
public interface IAtomKinetic extends IAtom {

    /**
     * Returns the velocity of the IAtom.  Modifying the returned IVector will
     * alter the IAtom's velocity.
     */
    public IVectorMutable getVelocity();
 
}
