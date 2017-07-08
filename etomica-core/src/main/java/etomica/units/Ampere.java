/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.units;
import java.io.ObjectStreamException;

import etomica.units.dimensions.Current;
import etomica.util.Constants;

/**
 * The ampere unit of electrical current.
 */
public final class Ampere extends SimpleUnit {

  /**
   * Singleton instance of this unit.
   */
    public static final Ampere UNIT = new Ampere();
    
    private Ampere() {
        //should check this more carefully
        // sqrt[ 1/(4 Pi epsilon0) * (A/m)^2 * (g/kg) * (D/g) * (A/m) * (s/ps)^3 ]
        super(Current.DIMENSION,
                Math.sqrt(1/4/Math.PI/8.854e-12*1e20*1000*Constants.AVOGADRO*1e10*1e-36), //2.326e11; conversion from Coulombs to (amu-A^3/ps^2)^(1/2)
                "amperes", "A", Prefix.ALLOWED);   
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
