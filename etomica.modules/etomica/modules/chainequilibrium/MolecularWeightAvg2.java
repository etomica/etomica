package etomica.modules.chainequilibrium;

import etomica.data.Data;
import etomica.data.DataPipe;
import etomica.data.DataProcessor;
import etomica.data.IDataInfo;
import etomica.data.types.DataDouble;
import etomica.data.types.DataDouble.DataInfoDouble;
import etomica.units.Quantity;

/**
 * Takes output from MeterChainLength as input and outputs the weight-average
 * molecular weight.
 * 
 * @author Andrew Schultz
 */
public class MolecularWeightAvg2 extends DataProcessor {

    public MolecularWeightAvg2() {
        data = new DataDouble();
        dataInfo = new DataInfoDouble("Avg MW", Quantity.DIMENSION);
    }

    protected Data processData(Data inputData) {
        double sum = 0, sum2 = 0;
        for (int i=0; i<inputData.getLength(); i++) {
            double v = inputData.getValue(i);
            sum += v * (i+1);
            sum2 += v;
        }
        data.x = sum / sum2;
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
