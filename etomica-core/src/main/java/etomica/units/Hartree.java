/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.units;
import etomica.units.dimensions.Energy;
import etomica.util.Constants;

public final class Hartree extends SimpleUnit {

    /**
     * Hartree unit of energy, corresponding to the approximate energy of
     * hydrogen in its ground state.
     */
    public static final Hartree UNIT = new Hartree();
  
    private Hartree() {
        super(Energy.DIMENSION,
        	2*Constants.PLANCK_H*Constants.RYDBERG_R*Constants.LIGHT_SPEED,  //262549.9617098284
        	"hartrees", "Ha", Prefix.ALLOWED
        	);
    }
}
