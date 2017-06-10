/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.potential;

import etomica.space.Space;

/**
 * Potential acting on 2 atoms or molecules.
 * 
 * @author David Kofke
 */

public abstract class Potential2 extends Potential {

    /**
     * Constructs potential with given space.
     */
    public Potential2(Space space) {
        super(2, space);
    }

}
