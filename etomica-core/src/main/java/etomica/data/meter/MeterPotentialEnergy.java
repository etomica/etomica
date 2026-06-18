/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.data.meter;

import etomica.box.Box;
import etomica.data.DataSourceScalar;
import etomica.data.IDataSourcePotential;
import etomica.potential.compute.PotentialCompute;
import etomica.units.dimensions.Energy;

public class MeterPotentialEnergy extends DataSourceScalar implements IDataSourcePotential {

    protected Box box;
    protected final PotentialCompute potentialMaster;
    protected boolean callComputeAll;

    public MeterPotentialEnergy(PotentialCompute potentialMaster) {
        super("Potential Energy", Energy.DIMENSION);
        this.potentialMaster = potentialMaster;
        callComputeAll = true;
    }

    /**
     * Computes total potential energy for box.
     * Currently, does not include long-range correction to truncation of energy
     */
    public double getDataAsScalar() {
        if (callComputeAll) {
            potentialMaster.computeAll(false);
        }
        return potentialMaster.getLastEnergy();
    }

    @Override
    public void doCallComputeAll(boolean callComputeAll) {
        this.callComputeAll = callComputeAll;
    }
}
