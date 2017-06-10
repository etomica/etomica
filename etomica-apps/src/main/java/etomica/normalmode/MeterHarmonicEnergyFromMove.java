/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.normalmode;

import etomica.data.DataSourceScalar;
import etomica.units.dimensions.Energy;

/**
 * Meter that returns the harmonic energy from the MC move as a data source.
 * 
 * @author Andrew Schultz
 */
public class MeterHarmonicEnergyFromMove extends DataSourceScalar {
    
    private static final long serialVersionUID = 1L;

    public MeterHarmonicEnergyFromMove(MCMoveHarmonic mcMoveHarmonic) {
        super("energy", Energy.DIMENSION);
        this.mcMoveHarmonic = mcMoveHarmonic;
    }

    public double getDataAsScalar() {
        return mcMoveHarmonic.getLastTotalEnergy();
    }

    protected final MCMoveHarmonic mcMoveHarmonic;
}
