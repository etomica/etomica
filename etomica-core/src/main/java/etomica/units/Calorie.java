/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.units;

import java.io.ObjectStreamException;

import etomica.units.dimensions.Energy;
import etomica.util.Constants;

/**
 * The Joule unit of energy, equal to 1 N-m or 1 kg-m^2/s^2.
 */
public final class Calorie extends SimpleUnit {

  /**
   * Singleton instance of this unit.
   */
    public static final Calorie UNIT = new Calorie();
    
    private Calorie() {
        super(Energy.DIMENSION,
        	Constants.AVOGADRO*1000.*1e20*1e-24*4.184, //6.022e22; conversion from calories to Joules kg-m^2/s^2 to Dalton-A^2/ps^2
        	"calories", "cal", Prefix.ALLOWED
        	);   
    }
}
