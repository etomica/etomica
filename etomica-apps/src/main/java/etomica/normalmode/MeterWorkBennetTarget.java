/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.normalmode;

import etomica.data.DataTag;
import etomica.data.IData;
import etomica.data.IDataSource;
import etomica.data.IDataInfo;
import etomica.data.meter.MeterPotentialEnergyFromIntegrator;
import etomica.data.types.DataDouble;
import etomica.data.types.DataDouble.DataInfoDouble;
import etomica.integrator.IntegratorBox;
import etomica.units.dimensions.Null;

/**
 * Meter used for direct sampling in the target-sampled system.  The meter
 * measures the energy difference between the Bennet's Overlap region and target
 * potentials.
 * 
 * @author Andrew Schultz & Tai Boon Tan
 */
public class MeterWorkBennetTarget implements IDataSource {
    
    public MeterWorkBennetTarget(IntegratorBox integrator, MeterHarmonicEnergy meterHarmonic) {
        meterTarget = new MeterPotentialEnergyFromIntegrator(integrator);
        this.integrator = integrator;
        this.meterHarmonic = meterHarmonic;
        
        data = new DataDouble();
        dataInfo = new DataInfoDouble("Target and Bennet's Overlap Energies", Null.DIMENSION);

        tag = new DataTag();
    }

    public IData getData() {
    	
    	double uTarget = meterTarget.getDataAsScalar();
    	double uHarmonic = meterHarmonic.getDataAsScalar();
    	double exp_uTarget = Math.exp(-(uTarget-latticeEnergy) /integrator.getTemperature());
    	double exp_uHarmonic = Math.exp(-uHarmonic /integrator.getTemperature());
    	double gamma = (exp_uTarget*exp_uHarmonic)/(exp_uTarget + refPref* exp_uHarmonic);
    	double overlapEnergy = -Math.log(gamma)*integrator.getTemperature();

//    	System.out.println("\nBennetTarget");
//    	System.out.println("uTarget-ulattice: "+ (uTarget - latticeEnergy));
//    	System.out.println("uHarmonic: "+ uHarmonic);
//    	System.out.println("uOverlap: "+overlapEnergy);
//    	System.out.println("uBennetTarget: "+((uTarget - latticeEnergy)-overlapEnergy));
        
    	data.x =  ((uTarget - latticeEnergy) - overlapEnergy)/integrator.getTemperature();
        return data;
    }

    public void setLatticeEnergy(double newLatticeEnergy) {
        latticeEnergy = newLatticeEnergy;
    }
    
    public IDataInfo getDataInfo() {
        return dataInfo;
    }

    public DataTag getTag() {
        return tag;
    }
	
    public double getRefPref() {
		return refPref;
	}

	public void setRefPref(double refPref) {
		this.refPref = refPref;
	}


    protected final MeterPotentialEnergyFromIntegrator meterTarget;
    protected final MeterHarmonicEnergy meterHarmonic;
    protected final IntegratorBox integrator;
    protected final DataDouble data;
    protected final DataInfoDouble dataInfo;
    protected final DataTag tag;
    protected double latticeEnergy;
    protected double refPref;

}
