package etomica.data.types;

import etomica.Data;
import etomica.DataInfo;
import etomica.data.DataGroupExtractor;
import etomica.data.DataJudge;
import etomica.data.DataProcessor;

/**
 * A DataProcessor that converts a Data instance into a DataArithmetic. Copies an
 * element of the input data to the DataDouble's encapsulated value and returns
 * the DataDouble instance.  Can cast from DataDouble, DataDoubleArray, DataInteger.
 * 
 * @author David Kofke
 *  
 */

/*
 * History Created on Jul 21, 2005 by kofke
 */
public class CastToArithmetic extends DataProcessor {

    protected DataInfo processDataInfo(DataInfo inputDataInfo) {
        Class inputClass = inputDataInfo.getClass();
        if (DataArithmetic.class.isAssignableFrom(inputClass)) {
            inputType = 0;
        } else if (inputClass == DataGroup.class) {
            inputType = 1;
            DataGroupExtractor extractor = new DataGroupExtractor(DataJudge.INSTANCEOF_ARITHMETIC);
            dataArithmetic = extractor.processData(null);
            
        } else {
            throw new IllegalArgumentException("Cannot cast to DataArithmetic from "
                    + inputClass);
        }
        return dataArithmetic.getDataInfo();
    }
    
    /**
     * Extracts a double from the input data and returns it encapsulated in a
     * DataDouble.
     * 
     * @param data
     *            a Data instance of the type indicated by the DataInfo at
     *            construction
     * @return a DataArithmetic holding the value cast from the given Data; the same
     *         instance is returned with every invocation.
     * 
     * @throws ClassCastException
     *             if the given Data is not of the same type as indicated by the
     *             DataInfo given at construction
     */
    protected Data processData(Data data) {
        switch (inputType) {
        case 0:
            return data;
        case 1:
            return dataArithmetic;
        default:
            throw new Error("Assertion error.  Input type out of range: "+inputType);
        }
    }
    
    public DataProcessor getDataCaster(DataInfo info) {
        return null;
    }

    private int inputType;
    private Data dataArithmetic;
}
