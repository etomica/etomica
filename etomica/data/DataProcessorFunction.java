package etomica.data;

import etomica.api.IData;
import etomica.api.IFunction;


/**
 * Applies a simple scalar function to all elements of the Data
 * that passes through.  Function is set at construction and
 * cannot be subsequently changed.
 *
 * @author David Kofke
 */
public class DataProcessorFunction extends DataProcessor {

    public DataProcessorFunction(IFunction function) {
        this.function = function;
    }
    
    /**
     * Applies the function to all elements of the input data.
     */
    protected IData processData(IData inputData) {
        data.E(inputData);
        data.map(function);
        return data;
    }

    /**
     * Returns a copy of the given dataInfo
     */
    protected IEtomicaDataInfo processDataInfo(IEtomicaDataInfo inputDataInfo) {
        dataInfo = inputDataInfo.getFactory().makeDataInfo();
        dataInfo.addTag(getTag());
        data = dataInfo.makeData();
        return dataInfo;
    }

    /**
     * Always returns null.
     */
    public DataPipe getDataCaster(IEtomicaDataInfo incomingDataInfo) {
        return null;
    }

    private static final long serialVersionUID = 1L;
    private final IFunction function;
    protected IData data;
}
