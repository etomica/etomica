/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.models.nitrogen;

import etomica.data.DataSourceScalar;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.data.meter.MeterPotentialEnergyFromIntegrator;
import etomica.integrator.IntegratorBox;
import etomica.units.dimensions.Null;

/**
 * Meter used for overlap sampling in the target-sampled system.  The meter
 * measures the ratio of the Boltzmann factors for the reference and target
 * potentials.
 * 
 * @author Tai Boon Tan
 */
public class MeterBoltzmannDirect extends DataSourceScalar {
	
	public MeterBoltzmannDirect(IntegratorBox integrator, MeterPotentialEnergy meterPotentialEnergy) {
    	super("Scaled Unit",Null.DIMENSION);
         
    	meterEnergy = new MeterPotentialEnergyFromIntegrator(integrator);
        this.integrator = integrator;
        this.meterPotentialEnergy = meterPotentialEnergy;
    }
    
	public double getDataAsScalar() {
		double uSampled = meterEnergy.getDataAsScalar();
    	double uMeasured = meterPotentialEnergy.getDataAsScalar();
    	
    	if(Double.isNaN(uSampled) || Double.isNaN(uMeasured)){
    		throw new RuntimeException("<MeterBoltzmannDirect> energy is NaN!!!!!!!!!!!!");
    	}
    	 
    	//System.out.println(uSampled + " " + uMeasured + " "+ (uMeasured-uSampled) + " " + Math.exp(-(uMeasured - uSampled) / integrator.getTemperature()));
    	return Math.exp(-(uMeasured - uSampled) / integrator.getTemperature());
	}
	
	private static final long serialVersionUID = 1L;
    protected final MeterPotentialEnergyFromIntegrator meterEnergy;
    protected final MeterPotentialEnergy meterPotentialEnergy;
    protected final IntegratorBox integrator;

}
