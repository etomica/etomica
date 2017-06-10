/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.units;

import etomica.units.dimensions.Angle;

import java.io.ObjectStreamException;

/**
 * Standard degree angular unit, such that for example a right angle is 90 degrees.
 */
public final class Degree extends SimpleUnit {

  /**
   * Single instance of this unit.
   */
    public static final Degree UNIT = new Degree();
    
    private Degree() {
        super(Angle.DIMENSION,
        	Math.PI/180., //conversion from degrees to radians
        	"degrees","\u00B0", Prefix.ALLOWED //unicode for the degree symbol
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
