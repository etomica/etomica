/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.units;
import etomica.units.dimensions.ElectricPotential;

import java.io.ObjectStreamException;

/**
 * The Volt unit of electric potential.  One volt is equal to a Joule per Coulomb,
 * or Watt per Ampere, or 1 meter^2 kg sec^-3 amp^-1.  One volt is approximately 25.885 
 * simulation units.
 */
public final class Volt extends SimpleUnit {

  /**
   * Singleton instance of this unit.
   */
    public static final Volt UNIT = new Volt();
    
    private Volt() {
        super(ElectricPotential.DIMENSION,
                Joule.UNIT.toSim(1.0)/Coulomb.UNIT.toSim(1.0), //25.885, conversion to (Angstrom^2 D/ps^2)*(4 pi eps0); 
	        "volts","V", Prefix.ALLOWED
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
