/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.normalmode;

import etomica.data.DataTag;
import etomica.data.IData;
import etomica.data.IDataSource;
import etomica.data.IDataInfo;
import etomica.data.types.DataDouble;
import etomica.data.types.DataDouble.DataInfoDouble;
import etomica.integrator.IntegratorBox;
import etomica.units.dimensions.Null;

/**
 * Meter used for compute for the probability of the harmonic sampling 
 * 
 * < e0 / (e1 + alpha * e0)>umb
 * 
 * @author Andrew Schultz & Tai Boon Tan
 */
public class MeterSamplingHarmonic implements IDataSource {
    
    public MeterSamplingHarmonic(IntegratorBox integrator, MCMoveAtomCoupledUmbrella move) {
        this.mcMove = move;
        this.integrator = integrator;
    	
        data = new DataDouble();
        dataInfo = new DataInfoDouble("Scaled Harmonic and soft sphere Energies", Null.DIMENSION);

        tag = new DataTag();
    }

    public IData getData() {
    	
    	double uTarget = mcMove.getPotentialEnergy();
    	double uHarmonic = mcMove.getHarmonicEnergy();
    	double exp_uTarget = Math.exp(-uTarget /integrator.getTemperature());
    	double exp_uHarmonic = Math.exp(-uHarmonic /integrator.getTemperature());
    	
       	double gamma = exp_uTarget + refPref*exp_uHarmonic;   	
    	data.x = exp_uHarmonic / gamma;
        
        return data;
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
    
    protected final MCMoveAtomCoupledUmbrella mcMove;
    protected final IntegratorBox integrator;
    protected final DataDouble data;
    protected final DataInfoDouble dataInfo;
    protected final DataTag tag;
    protected double refPref;



}
