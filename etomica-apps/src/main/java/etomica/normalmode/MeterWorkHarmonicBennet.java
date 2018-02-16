/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.normalmode;

import etomica.data.*;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.data.types.DataDouble;
import etomica.data.types.DataDouble.DataInfoDouble;
import etomica.potential.PotentialMaster;
import etomica.units.dimensions.Null;

/**
 * Meter used for overlap sampling in the harmonic-sampled system.  The meter
 * measures the ratio of the Boltzmann factors for the harmonic and target
 * potentials.
 * 
 * @author Andrew Schultz & Tai Boon Tan
 */
public class MeterWorkHarmonicBennet implements IDataSource {
    
    public MeterWorkHarmonicBennet(MCMoveHarmonic mcMoveHarmonic, PotentialMaster potentialMaster, double ref) {
        this.mcMoveHarmonic = mcMoveHarmonic;
        this.refPref = ref;
        meterEnergy = new MeterPotentialEnergy(potentialMaster, mcMoveHarmonic.getBox());
        data = new DataDouble();
        numSum = 0;
        denomSum = 0;
        dataInfo = new DataInfoDouble("Scaled Harmonic and soft sphere Energies", Null.DIMENSION);

        tag = new DataTag();
    }

    public void setTemperature(double newTemperature) {
        temperature = newTemperature;
    }
    
    public double getTemperature() {
        return temperature;
    }
    
    public IData getData() {
    	double e0 = Math.exp( -mcMoveHarmonic.getLastTotalEnergy()/ temperature);
    	double e1 = Math.exp( -(meterEnergy.getDataAsScalar() - latticeEnergy)/ temperature);
    	
    	double ratio = e0/(1+refPref*(e0/e1));     //e0*e1/(e1+refPref*e0);
    	double overlapEnergy = -Math.log(ratio);
    	double harmonicEnergy = mcMoveHarmonic.getLastTotalEnergy()/temperature;
        data.x = overlapEnergy - harmonicEnergy;
   
        /*
         * to take care of the limit when e1 --> 0
         */
        if (e1 < 1e-100){
            numSum += 0.0;              

        } else {
        	numSum += data.x*(ratio/e0);
        	
        }
        
        denomSum += ratio/e0;
        return data;
    }
    
    public Double getDataReweighted(){
    	
    	return numSum/denomSum;
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

    protected double temperature;
    protected final MeterPotentialEnergy meterEnergy;
    protected final MCMoveHarmonic mcMoveHarmonic;
    protected final DataDouble data;
    protected final DataInfoDouble dataInfo;
    protected final DataTag tag;
    protected double latticeEnergy;
    protected double numSum, denomSum;
    public double refPref;

}
