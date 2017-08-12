/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.units;
import etomica.units.dimensions.Dipole;

import java.io.ObjectStreamException;

/**
 * The debye unit of electrical dipole moment, equal to 10^-18 statC-cm,
 * or 3.33564e-30 C-m.
 *
 */
public final class Debye extends SimpleUnit {

  /**
   * Singleton instance of this unit.
   */
    public static final Debye UNIT = new Debye();
    
    private Debye() {
        super(Dipole.DIMENSION,
                //convert value in Coulomb-meter to (sim charge)-Angstrom
                Coulomb.UNIT.toSim(3.33564e-30)*1e10,//77.60 //Math.sqrt(Constants.AVOGADRO*1e40*1e-24)*1e-18, //77.6; conversion from (g-cm^5/s^2)^(1/2) to 10^-18 * (Daltons-A^5/ps^2)^(1/2) (Debye)
        	        "debyes",
        	        "D",
                 Prefix.ALLOWED);   
    }
}
