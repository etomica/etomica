/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.data.meter;

import etomica.data.DataSourceScalar;
import etomica.integrator.IntegratorMD;
import etomica.units.dimensions.Energy;

/**
 * Acts as a DataSource to retrieve the energy from the integrator
 */
public class MeterKineticEnergyFromIntegrator extends DataSourceScalar {

    public MeterKineticEnergyFromIntegrator(IntegratorMD aIntegrator) {
        super("Potential Energy", Energy.DIMENSION);
        integrator = aIntegrator;
    }

    public IntegratorMD getIntegrator() {
        return integrator;
    }

    public double getDataAsScalar() {
        return integrator.getKineticEnergy();
    }

    private final IntegratorMD integrator;
}
