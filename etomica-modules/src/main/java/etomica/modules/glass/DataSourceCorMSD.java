/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
package etomica.modules.glass;

import etomica.integrator.IntegratorMD;

public class DataSourceCorMSD extends DataSourceBlockAvgCor implements DataSourceMSD.MSDSink {

    // DataSourceMSD (by default) will collect interval<3, but not every time.  Our superclass
    // doesn't understand that.
    protected int minInterval = 3;

    public DataSourceCorMSD(IntegratorMD integrator) {
        super(integrator);
    }

    public void putMSD(int interval, long step, double msd) {
        if (interval < minInterval) return;
        processData(interval, step, msd);
    }
}
