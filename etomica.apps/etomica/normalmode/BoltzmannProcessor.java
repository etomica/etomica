package etomica.normalmode;

import etomica.data.Data;
import etomica.data.DataPipe;
import etomica.data.DataProcessor;
import etomica.data.IDataInfo;
import etomica.data.IDataInfoFactory;
import etomica.units.Null;
import etomica.util.Function;

/**
 * DataProcessor that returns the Boltzmann factor of the incoming energy.
 * @author Andrew Schultz
 */
public class BoltzmannProcessor extends DataProcessor {
    
    public BoltzmannProcessor() {
    }
    
    public IDataInfo processDataInfo(IDataInfo incomingDataInfo) {
        data = incomingDataInfo.makeData();
        IDataInfoFactory factory = incomingDataInfo.getFactory();
        // we get energy in, spit out unitless
        factory.setDimension(Null.DIMENSION);
        dataInfo = factory.makeDataInfo();
        return dataInfo;
    }
    
    public void setEnergyBase(double newEnergyBase) {
        energyBase = newEnergyBase;
    }
    
    public Data processData(Data incomingData) {
        data.E(incomingData);
        data.PE(-energyBase);
        data.TE(-1/temperature);
        data.map(Function.Exp.INSTANCE);
        return data;
    }
    
    public DataPipe getDataCaster(IDataInfo incomingDataInfo) {
        return null;
    }
    
    public void setTemperature(double newTemperature) {
        temperature = newTemperature;
    }
    
    public double getTemperature() {
        return temperature;
    }
    
    private static final long serialVersionUID = 1L;
    private Data data;
    private double temperature;
    private double energyBase;
}
