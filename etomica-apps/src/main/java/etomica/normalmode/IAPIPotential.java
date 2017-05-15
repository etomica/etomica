/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.normalmode;

import etomica.box.Box;

/**
 * Preliminary interface for proposed IPotential interface
 */
public interface IAPIPotential {

    /**
     * Returns the potential energy of the given box.
     */
    public double calculateEnergy(Box box);
}
