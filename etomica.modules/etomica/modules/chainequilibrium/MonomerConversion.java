package etomica.modules.chainequilibrium;

import etomica.data.Data;
import etomica.data.DataPipe;
import etomica.data.DataProcessor;
import etomica.data.IDataInfo;
import etomica.data.types.DataDouble;
import etomica.data.types.DataDouble.DataInfoDouble;
import etomica.units.Fraction;

/**
 * Takes output from MeterChainLength as input and outputs the monomer
 * conversion.
 * 
 * @author Andrew Schultz
 */
public class MonomerConversion extends DataProcessor {

    public MonomerConversion() {
        data = new DataDouble();
        dataInfo = new DataInfoDouble("Monomer conversion", Fraction.DIMENSION);
    }

    protected Data processData(Data inputData) {
        data.x = 1.0 - inputData.getValue(0);
        return data;
    }

    protected IDataInfo processDataInfo(IDataInfo inputDataInfo) {
        return dataInfo;
    }

    public DataPipe getDataCaster(IDataInfo inputDataInfo) {
        return null;
    }

    protected final DataDouble data;
    protected final DataInfoDouble dataInfo;
}
