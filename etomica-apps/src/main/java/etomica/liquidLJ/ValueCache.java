/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.liquidLJ;

import etomica.integrator.Integrator;
import etomica.data.DataSourceScalar;

public class ValueCache {
    protected long lastStep = -1;
    protected double lastValue;
    protected final DataSourceScalar dss;
    protected final Integrator integrator;
    public ValueCache(DataSourceScalar dss, Integrator integrator) {
        this.dss = dss;
        this.integrator = integrator;
    }
    public double getValue() {
        if (integrator.getStepCount() != lastStep) {
            lastStep = integrator.getStepCount();
            lastValue = dss.getDataAsScalar();
        }
        return lastValue;
    }
}
