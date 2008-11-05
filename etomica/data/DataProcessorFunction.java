package etomica.data;

import etomica.api.IData;
import etomica.api.IFunction;


/**
 * Applies a simple scalar function to all elements of the Data
 * that passes through.  Function is set at construction and
 * cannot be subsequently changed.
 *
 * @author David Kofke
 *
 */

public class DataProcessorFunction extends DataProcessor {

    public DataProcessorFunction(IFunction function) {
        this.function = function;
    }
    
    /**
     * Applies the function to all elements of the input data.
     * Uses the map method of the DataArithmetic interface.  
     * 
     * @throws ClassCastException if the input data does not implement DataArithmetic
     */
    protected IData processData(IData inputData) {
        inputData.map(function);
        return inputData;
    }

    /**
     * Returns the given DataInfo unchanged.
     * 
     * @throws IllegalArgumentException if the input data class does not 
     * implement DataArithmetic
     */
    protected IEtomicaDataInfo processDataInfo(IEtomicaDataInfo inputDataInfo) {
        dataInfo = inputDataInfo.getFactory().makeDataInfo();
        dataInfo.addTag(getTag());
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
}
