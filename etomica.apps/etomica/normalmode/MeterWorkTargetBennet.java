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
 * Meter used for overlap sampling in the target-sampled system.  The meter
 * measures the ratio of the Boltzmann factors for the harmonic and target
 * potentials.
 * 
 * @author Andrew Schultz & Tai Boon Tan
 */
public class MeterWorkTargetBennet implements IEtomicaDataSource {
    
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
      	//System.out.println("uHarmonic: "+ meterHarmonicEnergy.getDataAsScalar()+" uTarget: "+meterEnergy.getDataAsScalar()
		//			+" ,ratio: "+ratio+ " ,overlapEnergy: "+overlapEnergy);
    	denomSum += e0/(e1+refPref*e0);                //ratio/e1;
    	numSum +=   data.x*(e0/(e1+refPref*e0));      //data.x*(ratio/e1);
    	
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
