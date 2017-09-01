/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.units;

import etomica.units.dimensions.Length;

import java.io.ObjectStreamException;

/**
 * The meter unit of length, equal to 10^10 angstroms.
 */
public final class Meter extends SimpleUnit {

  /**
   * Singleton instance of this unit.
   */
    public static final Meter UNIT = new Meter();
    
    private Meter() {
    	    super(Length.DIMENSION,
    	            1e+10, //conversion from meters to Angstroms
    	            "meters", "m", Prefix.ALLOWED
    	    );
    }
}
