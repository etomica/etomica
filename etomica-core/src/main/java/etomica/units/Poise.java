/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.units;
import java.io.ObjectStreamException;

import etomica.units.dimensions.Viscosity;
import etomica.util.Constants;

/**
 * The Poise unit of viscosity, equal to 1 gram/(cm-sec).
 * This is the standard unit of viscosity in the CGS unit system.
 * It is equal to 0.1 pascal-seconds or approximately 6022.1 simulation viscosity units.
 */
public final class Poise extends SimpleUnit {

  /**
   * Singleton instance of this unit.
   */
    public static final Poise UNIT = new Poise();
  
  /**
   * Conversion factor to/from simulation units
   */
    private Poise() {
        super(Viscosity.DIMENSION,
                1e-8*1e-12*Constants.AVOGADRO, //6022.1; conversion from g/(cm-sec) to D/(A-ps)
                "poise", "P", Prefix.ALLOWED
        	);   
    }
}
