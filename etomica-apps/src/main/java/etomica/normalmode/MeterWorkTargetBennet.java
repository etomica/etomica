/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.normalmode;

import etomica.data.DataTag;
import etomica.data.IData;
import etomica.data.IDataSource;
import etomica.data.IEtomicaDataInfo;
import etomica.data.meter.MeterPotentialEnergyFromIntegrator;
import etomica.data.types.DataDouble;
import etomica.data.types.DataDouble.DataInfoDouble;
import etomica.integrator.IntegratorBox;
import etomica.units.dimensions.Null;

/**
 * Meter used for overlap sampling in the target-sampled system.  The meter
 * measures the ratio of the Boltzmann factors for the harmonic and target
 * potentials.
 * 
 * @author Andrew Schultz & Tai Boon Tan
 */
public class MeterWorkTargetBennet implements IDataSource {
    
    public MeterWorkTargetBennet(IntegratorBox integrator, MeterHarmonicEnergy meterHarmonicEnergy, double ref) {
        meterEnergy = new MeterPotentialEnergyFromIntegrator(integrator);
        this.integrator = integrator;
        this.meterHarmonicEnergy = meterHarmonicEnergy;
        this.refPref = ref;
        data = new DataDouble();
        numSum = 0;
        denomSum = 0;
        dataInfo = new DataInfoDouble("Scaled Harmonic and hard sphere Energies", Null.DIMENSION);

        tag = new DataTag();
    }

    public IData getData() {
    	double e0 = Math.exp(-meterHarmonicEnergy.getDataAsScalar()/integrator.getTemperature());
    	double e1 = Math.exp(-(meterEnergy.getDataAsScalar()-latticeEnergy)/integrator.getTemperature());
    		
    	double ratio = e0/(1+refPref*(e0/e1)); //e1*e0/(e1+refPref*e0);
    	double overlapEnergy = - Math.log(ratio);
        data.x = overlapEnergy - (meterEnergy.getDataAsScalar()-latticeEnergy)/integrator.getTemperature();

       	numSum +=   data.x*(e0/(e1+refPref*e0));      //data.x*(ratio/e1);
        denomSum += 1/((e1/e0)+refPref);                //ratio/e1;
        
        return data;
    }

    public double getDataReweighted(){
    	
    	return numSum/denomSum;
    }
    
    public void setLatticeEnergy(double newLatticeEnergy) {
        latticeEnergy = newLatticeEnergy;
    }
    
    public IEtomicaDataInfo getDataInfo() {
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

    protected final MeterPotentialEnergyFromIntegrator meterEnergy;
    protected final MeterHarmonicEnergy meterHarmonicEnergy;
    protected final IntegratorBox integrator;
    protected final DataDouble data;
    protected final DataInfoDouble dataInfo;
    protected final DataTag tag;
    protected double latticeEnergy;
    protected double numSum, denomSum;
    public double refPref;

}
