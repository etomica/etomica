/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.units;
import java.io.ObjectStreamException;

import etomica.units.dimensions.Quantity;
import etomica.util.Constants;

/**
 * The mole unit of quantity, approximately equal to 6.022e23 simulation units.
 */
public final class Mole extends SimpleUnit {

  /**
   * Singleton instance of this unit.
   */
    public static final Mole UNIT = new Mole();
    
    private Mole() {
        super(Quantity.DIMENSION,
               Constants.AVOGADRO, //6.022e23; conversion from moles to count (number)
               "moles",
               "mol",
               Prefix.ALLOWED
        	);
    }
}
