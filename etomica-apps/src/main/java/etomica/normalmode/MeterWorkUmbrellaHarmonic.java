/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.normalmode;

import etomica.data.DataTag;
import etomica.data.IData;
import etomica.data.IDataSource;
import etomica.data.IEtomicaDataInfo;
import etomica.data.types.DataDouble;
import etomica.data.types.DataDouble.DataInfoDouble;
import etomica.integrator.IntegratorBox;
import etomica.units.dimensions.Null;

/**
 * Meter used for direct sampling in the harmonic-sampled system.  The meter
 * measures the energy difference between the Umbrella's sampling region and harmonic system
 * 
 * @author Andrew Schultz & Tai Boon Tan
 */
public class MeterWorkUmbrellaHarmonic implements IDataSource {
    
    public MeterWorkUmbrellaHarmonic(IntegratorBox integrator, MCMoveAtomCoupledUmbrella move) {
        this.mcMove = move;
        this.integrator = integrator;
    	
        data = new DataDouble();
        dataInfo = new DataInfoDouble("Scaled Harmonic and soft sphere Energies", Null.DIMENSION);
        denomSum = 0;
        numSum = 0;
        tag = new DataTag();
    }

    public IData getData() {
    	
    	double uTarget = mcMove.getPotentialEnergy();
    	double uHarmonic = mcMove.getHarmonicEnergy();
    	double exp_uTarget = Math.exp(-(uTarget) /integrator.getTemperature());
    	double exp_uHarmonic = Math.exp(-uHarmonic /integrator.getTemperature());
//    	double exp_2uTarget = exp_uTarget*exp_uTarget;
//    	double exp_2uHarmonic = exp_uHarmonic*exp_uHarmonic;
    	
//    	double gamma = Math.sqrt(exp_2uTarget + refPref*refPref*exp_2uHarmonic);
    	double gamma = exp_uTarget + refPref*exp_uHarmonic;
    	double umbrellaEnergy = -Math.log(gamma)*integrator.getTemperature();

    	data.x = (uHarmonic- umbrellaEnergy)/integrator.getTemperature();
        
    	denomSum += exp_uHarmonic/gamma;
    	numSum += data.x*(exp_uHarmonic/gamma);    	
    	return data;
    }

    public double getDataReweighted(){
    	return numSum/denomSum;
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
    
    protected final MCMoveAtomCoupledUmbrella mcMove;
    protected final IntegratorBox integrator;
    protected final DataDouble data;
    protected final DataInfoDouble dataInfo;
    protected final DataTag tag;
    protected double numSum, denomSum;
    
    protected double refPref;



}
