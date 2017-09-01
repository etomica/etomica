/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.units;
import java.io.ObjectStreamException;

import etomica.units.dimensions.Charge;
import etomica.util.Constants;

/**
 * The Coulomb unit of electrical charge.
 */
public final class Coulomb extends SimpleUnit {

    /**
     * Singleton instance of this unit.
     */
    public static final Coulomb UNIT = new Coulomb();

    private Coulomb() {
        //note that Constants.EPSILON_0 is defined in terms of electron charge units, so must convert coulomb to e
        // then we have (1 electron-unit / 1.6e-19 C) / sqrt[4 Pi eps0]  
        super(Charge.DIMENSION,
                1.0 / (1.60217653e-19 * Math.sqrt(4 * Math.PI * Constants.EPSILON_0)), //2.326e21; conversion from Coulombs to (amu-A^3/ps^2)^(1/2)
                "coulombs", "C", Prefix.ALLOWED);
    }
}