package etomica.normalmode;

import etomica.api.IData;
import etomica.data.DataTag;
import etomica.data.IEtomicaDataInfo;
import etomica.data.IEtomicaDataSource;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.data.types.DataDouble;
import etomica.data.types.DataDouble.DataInfoDouble;
import etomica.potential.PotentialMaster;
import etomica.units.Null;

/**
 * Meter used for overlap sampling in the harmonic-sampled system.  The meter
 * measures the energy difference for the harmonic and target
 * potentials.
 * 
 * @author Andrew Schultz & Tai Boon Tan
 */
public class MeterWorkHarmonicPhaseSpace implements IEtomicaDataSource {
    
    public MeterWorkHarmonicPhaseSpace(MCMoveHarmonic mcMoveHarmonic, PotentialMaster potentialMaster) {
        this.mcMoveHarmonic = mcMoveHarmonic;
        meterEnergy = new MeterPotentialEnergy(potentialMaster);
        meterEnergy.setBox(mcMoveHarmonic.getBox());
        data = new DataDouble();
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
        data.x = ((meterEnergy.getDataAsScalar() - latticeEnergy) -
        		mcMoveHarmonic.getLastTotalEnergy()) / temperature;
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

    protected double temperature;
    protected final MeterPotentialEnergy meterEnergy;
    protected final MCMoveHarmonic mcMoveHarmonic;
    protected final DataDouble data;
    protected final DataInfoDouble dataInfo;
    protected final DataTag tag;
    protected double latticeEnergy;
}
