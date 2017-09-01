/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.units;
import etomica.units.dimensions.Temperature;
import etomica.util.Constants;

public final class Kelvin extends SimpleUnit {

  /**
   * Singleton instance of this unit.
   */
    public static final Kelvin UNIT = new Kelvin();
  
    private Kelvin() {
        super(Temperature.DIMENSION,
        	Constants.BOLTZMANN_K,//convert to simulation energy units by multiplying by Boltzmann's constant
        	"kelvins", "K", Prefix.ALLOWED
        	);
    }
}
