/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.data.meter;

import etomica.data.DataSourceScalar;
import etomica.integrator.IntegratorBoxFasterer;
import etomica.units.dimensions.Energy;

/**
 * Acts as a DataSource to retrieve the energy from the integrator
 */
public class MeterPotentialEnergyFromIntegratorFasterer extends DataSourceScalar {

    public MeterPotentialEnergyFromIntegratorFasterer() {
        super("Potential Energy", Energy.DIMENSION);
    }

    public MeterPotentialEnergyFromIntegratorFasterer(IntegratorBoxFasterer aIntegrator) {
        this();
        setIntegrator(aIntegrator);
    }

    public void setIntegrator(IntegratorBoxFasterer newIntegrator) {
        integrator = newIntegrator;
    }

    public IntegratorBoxFasterer getIntegrator() {
        return integrator;
    }

    public double getDataAsScalar() {
        return integrator.getPotentialEnergy();
    }

    private IntegratorBoxFasterer integrator;
}
