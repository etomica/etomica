package etomica.normalmode;

import etomica.data.Data;
import etomica.data.DataPipe;
import etomica.data.DataProcessor;
import etomica.data.IDataInfo;
import etomica.data.types.DataDouble;

/**
 * DataProcessor that sends on the sum of the incoming data values as a scalar.
 * @author Andrew Schultz
 */
public class DataProcessorSum extends DataProcessor {
    public IDataInfo processDataInfo(IDataInfo incomingDataInfo) {
        data = new DataDouble();
        dataInfo = new DataDouble.DataInfoDouble(incomingDataInfo.getLabel(), incomingDataInfo.getDimension());
        return dataInfo;
    }
    
    public Data processData(Data incomingData) {
        data.x = 0;
        int n = incomingData.getLength();
        for (int i=0; i<n; i++) {
            data.x += incomingData.getValue(i);
        }
        return data;
    }
    
    public DataPipe getDataCaster(IDataInfo incomingDataInfo) {
        return null;
    }
    
    private static final long serialVersionUID = 1L;
    private DataDouble data;
}
