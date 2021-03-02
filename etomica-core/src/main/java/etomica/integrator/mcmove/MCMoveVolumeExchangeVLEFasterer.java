/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.integrator.mcmove;

import etomica.data.meter.MeterDensity;
import etomica.integrator.IntegratorBoxFasterer;
import etomica.space.Space;
import etomica.util.random.IRandom;

public class MCMoveVolumeExchangeVLEFasterer extends MCMoveVolumeExchangeFasterer {
    private final MeterDensity meterDensity1;
    private final MeterDensity meterDensity2;


    public MCMoveVolumeExchangeVLEFasterer(IRandom random, Space space,
                                           IntegratorBoxFasterer integrator1, IntegratorBoxFasterer integrator2) {
        super(random, space, integrator1, integrator2);
        meterDensity1 = new MeterDensity(integrator1.getBox());
        meterDensity2 = new MeterDensity(integrator2.getBox());
    }

    public boolean doTrial() {
        boolean success = super.doTrial();
        if (!success) return false;
        double density1 = meterDensity1.getDataAsScalar();
        double density2 = meterDensity2.getDataAsScalar();
        if (density2 > density1) {
            rejectNotify();
            return false;
        }
        return true;
    }
}
