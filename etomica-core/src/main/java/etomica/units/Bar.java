/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.units;
import java.io.ObjectStreamException;

import etomica.units.dimensions.Pressure;
import etomica.util.Constants;

/**
 * The bar unit of pressure, equal to 10^5 N/m^2.
 * Equal to approximately 0.006022 simulation pressure units.
 */
public final class Bar extends SimpleUnit {

  /**
   * Singleton instance of this unit.
   */
    public static final Bar UNIT = new Bar();
  
  /**
   * Conversion factor to/from simulation units
   */
    private Bar() {
        super(Pressure.DIMENSION,
                1e5*1000.*1e-10*1e-24*Constants.AVOGADRO, //6.022e-3; conversion from 10^5 kg/(m-s^2) to amu/(A-ps^2)
                "bars", "bar", Prefix.ALLOWED
        	);   
    }
}
