/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.data.meter;

import etomica.data.DataSourceScalar;
import etomica.integrator.IntegratorBox;
import etomica.units.dimensions.Energy;

/**
 * Acts as a DataSource to retrieve the energy from the integrator 
 */
public class MeterPotentialEnergyFromIntegrator extends DataSourceScalar {

    public MeterPotentialEnergyFromIntegrator() {
        super("Potential Energy",Energy.DIMENSION);
    }
    
    public MeterPotentialEnergyFromIntegrator(IntegratorBox aIntegrator) {
        this();
        setIntegrator(aIntegrator);
    }
    
    public void setIntegrator(IntegratorBox newIntegrator) {
        integrator = newIntegrator;
    }
    
    public IntegratorBox getIntegrator() {
        return integrator;
    }
    
    public double getDataAsScalar() {
        return integrator.getPotentialEnergy();
    }
    
    private static final long serialVersionUID = 1L;
    private IntegratorBox integrator;
}
