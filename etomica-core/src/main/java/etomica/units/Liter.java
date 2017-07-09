/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.units;

import etomica.units.dimensions.Volume;

import java.io.ObjectStreamException;

/**
 * The liter unit of volume, equal to 1000 cm^3 or 0.001 m^3 or 10^27 A^3.
 */
public final class Liter extends SimpleUnit {

  /**
   * Singleton instance of this class.
   */
    public static final Liter UNIT = new Liter();
    
    public Liter() {
        super(Volume.DIMENSION,
        	1e+27, //conversion from liters to Angstroms^3 (10^-3 m^3/liter * (10^10 Angstroms/meter)^3)
        	"liters", "l", Prefix.ALLOWED
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
