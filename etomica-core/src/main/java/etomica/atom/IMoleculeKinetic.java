/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.atom;

import etomica.space.Vector;


/**
 * Interface for an molecule that holds vectors for velocity.
 */
public interface IMoleculeKinetic extends IMolecule {

    public Vector getVelocity();
 
}
