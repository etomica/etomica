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

    public MeterKineticEnergyFromIntegrator() {
        super("Kinetic Energy",Energy.DIMENSION);
    }
    
    public MeterKineticEnergyFromIntegrator(IntegratorMD aIntegrator) {
        this();
        setIntegrator(aIntegrator);
    }
    
    public void setIntegrator(IntegratorMD newIntegrator) {
        integrator = newIntegrator;
    }
    
    public IntegratorMD getIntegrator() {
        return integrator;
    }
    
    public double getDataAsScalar() {
        return integrator.getKineticEnergy();
    }
 
    private static final long serialVersionUID = 1L;
    private IntegratorMD integrator;
}
