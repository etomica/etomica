/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.units;
import etomica.units.dimensions.Pressure;

import java.io.ObjectStreamException;

/**
 * The barye unit of pressure, equal to 1 dyn/cm^2.
 * This is the standard unit of pressure in the CGS unit system.
 * It is equal to 0.1 bar or 6.022e-4 simulation pressure units.
 */
public final class Barye extends SimpleUnit {

  /**
   * Singleton instance of this unit.
   */
    public static final Barye UNIT = new Barye();
  
  /**
   * Conversion factor to/from simulation units
   */
    private Barye() {
        super(Pressure.DIMENSION,
                Bar.UNIT.toSim(0.1), //6.022e-4; conversion from 10^5 kg/(m-s^2) to amu/(A-ps^2)
                "baryes", "Ba", Prefix.ALLOWED
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
