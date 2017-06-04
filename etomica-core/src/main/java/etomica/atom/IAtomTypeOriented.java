/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.atom;

import etomica.space.Vector;

public interface IAtomTypeOriented extends IAtomType {

    /**
     * Returns the principle components of the moment of inertia of the
     * atom within the body-fixed frame.  Do NOT modify the returned moment
     * of inertia returned.
     */
    public Vector getMomentOfInertia();
}
