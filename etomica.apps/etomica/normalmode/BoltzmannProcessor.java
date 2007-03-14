package etomica.normalmode;

import etomica.data.Data;
import etomica.data.DataInfo;
import etomica.data.DataInfoFactory;
import etomica.data.DataProcessor;
import etomica.units.Null;
import etomica.util.Function;

/**
 * DataProcessor that returns the Boltzmann factor of the incoming energy.
 * @author Andrew Schultz
 */
public class BoltzmannProcessor extends DataProcessor {
    public DataInfo processDataInfo(DataInfo incomingDataInfo) {
        data = incomingDataInfo.makeData();
        DataInfoFactory factory = incomingDataInfo.getFactory();
        // we get energy in, spit out unitless
        factory.setDimension(Null.DIMENSION);
        dataInfo = factory.makeDataInfo();
        return dataInfo;
    }
    
    public Data processData(Data incomingData) {
        data.E(incomingData);
        data.TE(-1/temperature);
        data.map(Function.Exp.INSTANCE);
        return data;
    }
    
    public DataProcessor getDataCaster(DataInfo incomingDataInfo) {
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
}
