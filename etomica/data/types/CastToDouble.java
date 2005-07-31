package etomica.data.types;

import etomica.Data;
import etomica.DataInfo;
import etomica.data.DataGroupExtractor;
import etomica.data.DataJudge;
import etomica.data.DataProcessor;

/**
 * A DataProcessor that converts a Data instance into a DataDouble. Copies an
 * element of the input data to the DataDouble's encapsulated value and returns
 * the DataDouble instance.  Can cast from DataDouble, DataDoubleArray, DataInteger,
 * DataGroup.
 * 
 * @author David Kofke
 *  
 */

/*
 * History Created on Jul 21, 2005 by kofke
 */
public class CastToDouble extends DataProcessor {

    protected DataInfo processDataInfo(DataInfo inputDataInfo) {
        dataDouble = new DataDouble(inputDataInfo.getLabel(), inputDataInfo.getDimension());
        Class inputClass = inputDataInfo.getDataClass();
        if (inputClass == DataDouble.class) {
            inputType = 0;
        } else if (inputClass == DataDoubleArray.class) {
            inputType = 1;
        } else if (inputClass == DataInteger.class) {
            inputType = 2;
        } else if (inputClass == DataGroup.class) {
            inputType = 3;
            DataGroupExtractor extractor = new DataGroupExtractor(new DataJudge.ByClass(DataDouble.class, true));
            dataDouble = (DataDouble)extractor.processData(null);
        } else {
            throw new IllegalArgumentException("Cannot cast to double from "
                    + inputClass);
        }
        return dataDouble.getDataInfo();
    }
    
    /**
     * Extracts a double from the input data and returns it encapsulated in a
     * DataDouble.
     * 
     * @param data
     *            a Data instance of the type indicated by the DataInfo at
     *            construction
     * @return a DataDouble holding the value cast from the given Data; the same
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
            dataDouble.x = ((DataDoubleArray) data).getData()[0];
            return dataDouble;
        case 2:
            dataDouble.x = ((DataInteger) data).x;
            return dataDouble;
        case 3:
            return dataDouble;
        default:
            throw new Error("Assertion error.  Input type out of range: "+inputType);
        }
    }
    
    public DataProcessor getDataCaster(DataInfo info) {
        return null;
    }

    private DataDouble dataDouble;
    private int inputType;
}
