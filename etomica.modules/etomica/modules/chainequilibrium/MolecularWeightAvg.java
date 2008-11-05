package etomica.modules.chainequilibrium;

import etomica.api.IData;
import etomica.data.DataPipe;
import etomica.data.DataProcessor;
import etomica.data.IEtomicaDataInfo;
import etomica.data.types.DataDouble;
import etomica.data.types.DataDouble.DataInfoDouble;
import etomica.units.Quantity;

/**
 * Takes output from MeterChainLength as input and outputs the number-average
 * molecular weight.
 * 
 * @author Andrew Schultz
 */
public class MolecularWeightAvg extends DataProcessor {

    public MolecularWeightAvg() {
        data = new DataDouble();
        dataInfo = new DataInfoDouble("Avg MW", Quantity.DIMENSION);
    }

    protected IData processData(IData inputData) {
        double sum = 0, sum2 = 0;
        for (int i=0; i<inputData.getLength(); i++) {
            double v = inputData.getValue(i);
            sum += v;
            sum2 += v / (i+1);
        }
        data.x = sum / sum2;
        return data;
    }

    protected IEtomicaDataInfo processDataInfo(IEtomicaDataInfo inputDataInfo) {
        return dataInfo;
    }

    public DataPipe getDataCaster(IEtomicaDataInfo inputDataInfo) {
        return null;
    }

    protected final DataDouble data;
    protected final DataInfoDouble dataInfo;
}
