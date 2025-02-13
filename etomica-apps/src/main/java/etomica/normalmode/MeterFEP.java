/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.normalmode;

import etomica.box.Box;
import etomica.data.DataSourceScalar;
import etomica.data.IDataSourcePotential;
import etomica.potential.compute.PotentialCompute;
import etomica.units.dimensions.Energy;

public class MeterFEP extends DataSourceScalar implements IDataSourcePotential {

    protected Box box;
    protected final PotentialCompute potentialMaster;
    protected boolean callComputeAll;
    protected double U0, beta;
    public MeterFEP(PotentialCompute potentialMaster, double beta) {
        super("Potential Energy", Energy.DIMENSION);
        this.potentialMaster = potentialMaster;
        callComputeAll = true;
        this.U0 = 0;
        this.beta = beta;
    }

    public void setU0(double U0) {
        this.U0 = U0;
    }

    /**
     * Computes total potential energy for box.
     * Currently, does not include long-range correction to truncation of energy
     */
    public double getDataAsScalar() {
        if (callComputeAll) {
            potentialMaster.computeAll(false);
        }
        double dU = potentialMaster.getLastEnergy() - U0;
        double expU = Math.exp(beta*dU);
        return expU;
    }

    @Override
    public void doCallComputeAll(boolean callComputeAll) {
        this.callComputeAll = callComputeAll;
    }
}
