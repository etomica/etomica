package etomica.normalmode;

import etomica.api.IData;
import etomica.data.DataTag;
import etomica.data.IEtomicaDataInfo;
import etomica.data.IEtomicaDataSource;
import etomica.data.meter.MeterPotentialEnergyFromIntegrator;
import etomica.data.types.DataDouble;
import etomica.data.types.DataDouble.DataInfoDouble;
import etomica.integrator.IntegratorBox;
import etomica.units.Null;

/**
 *Meter used for compute for the probability of the target system sampling 
 * 
 * < e1 / sqrt(e1^2 + alpha^2 * e0^2)>1
 * 
 * @author Andrew Schultz & Tai Boon Tan
 */
public class MeterSamplingTarget implements IEtomicaDataSource {
    
    public MeterSamplingTarget(IntegratorBox integrator, MeterHarmonicEnergy meterHarmonic) {
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
    	double exp_2uTarget = exp_uTarget*exp_uTarget;
    	double exp_2uHarmonic = exp_uHarmonic*exp_uHarmonic;
    	
    	double gamma = Math.sqrt(exp_2uTarget + refPref*refPref*exp_2uHarmonic);

//    	System.out.println("\nBennetTarget");
//    	System.out.println("uTarget-ulattice: "+ (uTarget - latticeEnergy));
//    	System.out.println("uHarmonic: "+ uHarmonic);
//    	System.out.println("uOverlap: "+overlapEnergy);
//    	System.out.println("uBennetTarget: "+((uTarget - latticeEnergy)-overlapEnergy));
        
    	data.x =  exp_uTarget /gamma;
        return data;
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


    protected final MeterPotentialEnergyFromIntegrator meterTarget;
    protected final MeterHarmonicEnergy meterHarmonic;
    protected final IntegratorBox integrator;
    protected final DataDouble data;
    protected final DataInfoDouble dataInfo;
    protected final DataTag tag;
    protected double latticeEnergy;
    protected double refPref;

}
