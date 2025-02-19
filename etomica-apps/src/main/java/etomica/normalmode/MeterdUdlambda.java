/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.normalmode;

import etomica.box.Box;
import etomica.data.DataSourceScalar;
import etomica.data.IDataSourcePotential;
import etomica.potential.compute.PotentialCompute;
import etomica.units.dimensions.Energy;

public class MeterdUdlambda extends DataSourceScalar implements IDataSourcePotential {

    protected Box box;
    protected final PotentialCompute pm1, pm2;
    protected boolean callComputeAll;
    public MeterdUdlambda(PotentialCompute pm1, PotentialCompute pm2) {
        super("dU/dlambda", Energy.DIMENSION);
        this.pm1 = pm1;
        this.pm2 = pm2;
        callComputeAll = true;
    }

    public double getDataAsScalar() {
        if (callComputeAll) {
            pm1.computeAll(false);
            pm2.computeAll(false);
        }
        double dUdl = pm1.getLastEnergy() - pm2.getLastEnergy();
        return dUdl;
    }

    @Override
    public void doCallComputeAll(boolean callComputeAll) {
        this.callComputeAll = callComputeAll;
    }
}
