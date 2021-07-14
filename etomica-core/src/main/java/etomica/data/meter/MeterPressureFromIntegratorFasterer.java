/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.data.meter;

import etomica.box.Box;
import etomica.data.DataSourceScalar;
import etomica.integrator.IntegratorBoxFasterer;
import etomica.units.dimensions.Pressure;
import etomica.units.dimensions.Pressure2D;

/**
 * Acts as a DataSource to retrieve the pressure from the integrator's potential master
 */
public class MeterPressureFromIntegratorFasterer extends DataSourceScalar {

    public MeterPressureFromIntegratorFasterer(IntegratorBoxFasterer aIntegrator) {
        super("Pressure", aIntegrator.getBox().getSpace().D() == 3 ? Pressure.DIMENSION : Pressure2D.DIMENSION);
        integrator = aIntegrator;
    }

    public IntegratorBoxFasterer getIntegrator() {
        return integrator;
    }

    public double getDataAsScalar() {
        double lastVirial = integrator.getPotentialCompute().getLastVirial();
        Box box = integrator.getBox();
        double vol = box.getBoundary().volume();
        int dim = box.getSpace().D();
        double temperature = integrator.getTemperature();
        return (box.getMoleculeList().size() / vol) * temperature
                - lastVirial / (vol * dim);
    }

    private final IntegratorBoxFasterer integrator;
}
