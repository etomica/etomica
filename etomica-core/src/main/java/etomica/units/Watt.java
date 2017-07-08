/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.units;

import etomica.units.dimensions.Energy;
import etomica.units.dimensions.Power;
import etomica.util.Constants;

import java.io.ObjectStreamException;

/**
 * The Watt unit of power, equal to 1 J/s or 1 kg-m^2/s^3.
 */
public final class Watt extends SimpleUnit {

  /**
   * Singleton instance of this unit.
   */
    public static final Watt UNIT = new Watt();

    private Watt() {
        super(Power.DIMENSION,
        	Constants.AVOGADRO*1000.*1e20*1e-36, //6.022e22; conversion from kg-m^2/s^3 to Dalton-A^2/ps^3
        	"watts", "W", Prefix.ALLOWED
        	);   
    }
}
