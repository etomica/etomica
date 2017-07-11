/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.units;
import java.io.ObjectStreamException;

import etomica.units.dimensions.Energy;
import etomica.util.Constants;

/**
 * The electronvolt unit of energy, equal to approximately 1.602e-19 Joules.
 */
public final class ElectronVolt extends SimpleUnit {

  /**
   * Singleton instance of this unit.
   */
    public static final ElectronVolt UNIT = new ElectronVolt();

    private ElectronVolt() {
        super(Energy.DIMENSION,
         // We first convert the eV into J, and then the J into Da*Angstrom^2/(ps^2)
         // (1.60217653e-19 J/eV)*(1 kg-m^2/s^2/J) (1000 g/kg) (N_avo D/g) (1e10 A/m)^2 (1e-12 s/ps)^2
        	1.60217653e-19 * Constants.AVOGADRO * 1000. * 1e20 * 1e-24, // conversion from eV to Dalton-A^2/ps^2
        	"ElectronVolts", "eV", Prefix.ALLOWED
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
