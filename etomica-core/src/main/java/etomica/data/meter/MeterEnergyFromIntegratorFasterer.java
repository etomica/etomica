/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.data.meter;

import etomica.data.DataSourceScalar;
import etomica.integrator.IntegratorMDFasterer;
import etomica.units.dimensions.Energy;

/**
 * Acts as a DataSource to retrieve the energy from the integrator
 */
public class MeterEnergyFromIntegratorFasterer extends DataSourceScalar {

    public MeterEnergyFromIntegratorFasterer(IntegratorMDFasterer aIntegrator) {
        super("Total Energy", Energy.DIMENSION);
        integrator = aIntegrator;
    }

    public IntegratorMDFasterer getIntegrator() {
        return integrator;
    }

    public double getDataAsScalar() {
        return integrator.getKineticEnergy() + integrator.getPotentialEnergy();
    }

    private final IntegratorMDFasterer integrator;
}
