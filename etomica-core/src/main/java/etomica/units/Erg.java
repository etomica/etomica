/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.units;
import etomica.units.dimensions.Energy;

import java.io.ObjectStreamException;

/**
 * The erg unit of energy, equal to 1 dyn-cm or 1 g-cm^2/s^2.
 * This is the standard unit of energy in the CGS unit system.
 * It is equal to 10^-7 Joules, or approximately 6.022e15 simulation energy units.
 */
public final class Erg extends SimpleUnit {

  /**
   * Singleton instance of this unit.
   */
    public static final Erg UNIT = new Erg();
    
    private Erg() {
        super(Energy.DIMENSION,
        	Joule.UNIT.toSim(1.e-7), //6.022e15; conversion from g-cm^2/s^2 to Dalton-A^2/ps^2
        	"ergs", "erg", Prefix.ALLOWED
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
