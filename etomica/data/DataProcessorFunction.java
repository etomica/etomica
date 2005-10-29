package etomica.data;

import etomica.data.types.DataArithmetic;
import etomica.util.Function;


/**
 * Applies a simple scalar function to all elements of the Data
 * that passes through.  Function is set at construction and
 * cannot be subsequently changed.
 *
 * @author David Kofke
 *
 */

/*
 * History
 * Created on Oct 29, 2005 by kofke
 */
public class DataProcessorFunction extends DataProcessor {

    public DataProcessorFunction(Function function) {
        this.function = function;
    }
    
    /**
     * Applies the function to all elements of the input data.
     * Uses the map method of the DataArithmetic interface.  
     * 
     * @throws ClassCastException if the input data does not implement DataArithmetic
     */
    protected Data processData(Data inputData) {
        ((DataArithmetic)inputData).map(function);
        return inputData;
    }

    /**
     * Returns the given DataInfo unchanged.
     */
    protected DataInfo processDataInfo(DataInfo inputDataInfo) {
        return inputDataInfo;
    }

    /**
     * Always returns null.
     */
    public DataProcessor getDataCaster(DataInfo dataInfo) {
        return null;
    }

    private final Function function;
}
