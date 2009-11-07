package etomica.normalmode;

import etomica.data.DataTag;
import etomica.data.IData;
import etomica.data.IEtomicaDataInfo;
import etomica.data.IEtomicaDataSource;
import etomica.data.types.DataDouble;
import etomica.data.types.DataDouble.DataInfoDouble;
import etomica.integrator.IntegratorBox;
import etomica.units.Null;

/**
 *Meter used for compute for the probability of the target system sampling 
 * 
 * < e1 / (e1 + alpha * e0)>umb
 * 
 * @author Tai Boon Tan
 */
public class MeterSamplingTarget implements IEtomicaDataSource {
    
    public MeterSamplingTarget(IntegratorBox integrator, MCMoveAtomCoupledUmbrella move) {
        this.mcMove = move;
        this.integrator = integrator;
        
        data = new DataDouble();
        dataInfo = new DataInfoDouble("Target and Umbrella Energies", Null.DIMENSION);

        tag = new DataTag();
    }

    public IData getData() {
    	
    	double uTarget = mcMove.getPotentialEnergy();
    	double uHarmonic = mcMove.getHarmonicEnergy();
    	double exp_uTarget = Math.exp(-uTarget /integrator.getTemperature());
    	double exp_uHarmonic = Math.exp(-uHarmonic /integrator.getTemperature());
//    	double exp_2uTarget = exp_uTarget*exp_uTarget;
//    	double exp_2uHarmonic = exp_uHarmonic*exp_uHarmonic;
    	
//    	double gamma = Math.sqrt(exp_2uTarget + refPref*refPref*exp_2uHarmonic);
    	double gamma = exp_uTarget + refPref*exp_uHarmonic;

        
    	data.x =  exp_uTarget /gamma;
        return data;
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
    protected double refPref;

}
