/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.units;
import java.io.ObjectStreamException;

import etomica.units.dimensions.Pressure;
import etomica.util.Constants;

/**
 * The Pascal unit of pressure, equal to 1 N/m^2.
 */
public final class Pascal extends SimpleUnit {

  /**
   * Singleton instance of this unit.
   */
    public static final Pascal UNIT = new Pascal();
  
  /**
   * Conversion factor to/from simulation units
   */
    private Pascal() {
        super(Pressure.DIMENSION,
                Constants.AVOGADRO*1000.*1e-10*1e-24, //6.022e-8; conversion from kg/(m-s^2) to D/(A-ps^2)
                "pascals", "Pa", Prefix.ALLOWED
        	);   
    }
}
