package etomica.normalmode;

import etomica.data.Data;
import etomica.data.DataSource;
import etomica.data.DataTag;
import etomica.data.IDataInfo;
import etomica.data.meter.MeterPotentialEnergyFromIntegrator;
import etomica.data.types.DataDouble;
import etomica.data.types.DataDouble.DataInfoDouble;
import etomica.integrator.IntegratorBox;
import etomica.units.Null;

/**
 * Meter used for overlap sampling in the target-sampled system.  The meter
 * measures the boltzmann factor difference for the harmonic and target
 * potentials.
 * 
 * @author Andrew Schultz
 */
public class MeterDirectSamplingTarget implements DataSource {
    
    public MeterDirectSamplingTarget(IntegratorBox integrator, MeterHarmonicEnergy meterHarmonicEnergy) {
        meterEnergy = new MeterPotentialEnergyFromIntegrator(integrator);
        this.integrator = integrator;
        this.meterHarmonicEnergy = meterHarmonicEnergy;
        data = new DataDouble();
        dataInfo = new DataInfoDouble("Scaled Harmonic and soft sphere Energies", Null.DIMENSION);

        tag = new DataTag();
    }

    public Data getData() {
        data.x = Math.exp(-(meterHarmonicEnergy.getDataAsScalar() -
                (meterEnergy.getDataAsScalar() - latticeEnergy)) / integrator.getTemperature());
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

    protected final MeterPotentialEnergyFromIntegrator meterEnergy;
    protected final MeterHarmonicEnergy meterHarmonicEnergy;
    protected final IntegratorBox integrator;
    protected final DataDouble data;
    protected final DataInfoDouble dataInfo;
    protected final DataTag tag;
    protected double latticeEnergy;
}
