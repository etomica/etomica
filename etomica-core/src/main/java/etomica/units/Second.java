/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.units;

import etomica.units.dimensions.Time;

import java.io.ObjectStreamException;

public final class Second extends SimpleUnit {

  /**
   * Singleton instance of this unit.
   */
    public static final Second UNIT = new Second();
    
    private Second() {
        super(Time.DIMENSION,
                1e+12, //conversion from seconds to picoseconds
                "seconds", "s", Prefix.ALLOWED
        	);   
    }
    
    /**
     * Required to guarantee singleton when deserializing.
     * 
     * @return the singleton UNIT
     */
    private Object readResolve() throws ObjectStreamException {
        return UNIT;
    }
    
    private static final long serialVersionUID = 1;

}
