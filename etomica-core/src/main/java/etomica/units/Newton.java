/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.units;
import java.io.ObjectStreamException;

import etomica.units.dimensions.Force;
import etomica.util.Constants;

/**
 * The Newton unit of force, equal to 1 kg-m/s^2.
 */
public final class Newton extends SimpleUnit {

  /**
   * Singleton instance of this unit.
   */
    public static final Newton UNIT = new Newton();
    
    private Newton() {
        super(Force.DIMENSION,
        	Constants.AVOGADRO*(1000.*1e10*1e-24), //6.022e12; conversion from kg-m/s^2 to Dalton-A/ps^2
        	"newtons", "N", Prefix.ALLOWED
        	);   
    }
}
