package etomica.data.types;

import etomica.data.Data;
import etomica.data.DataInfo;
import etomica.data.DataProcessor;
import etomica.data.DataTag;
import etomica.data.types.DataDouble.DataInfoDouble;
import etomica.data.types.DataDoubleArray.DataInfoDoubleArray;
import etomica.data.types.DataInteger.DataInfoInteger;

/**
 * A DataProcessor that converts a Data instance into a DataInteger. Copies an
 * element of the input data to the DataIntegers's encapsulated value and returns
 * the DataInteger instance.
 * Casting for various types of Data is performed as follows:
 * <ul>
 * <li><u>DataDoubleArray</u>. Uses only first element of array, cast to <tt>int</tt>.
 * 
 * <li><u>DataDouble</u>. Casts the encapsulated <tt>double</tt> to an <tt>int</tt>
 * 
 * <li><u>DataInteger</u>. Does nothing, passing the Data through.
 * 
 * <li><u>DataGroup</u>. Uses DataExtractor to locate a DataInteger in DataGroup.  Other 
 * Data in group are discarded. 
 * </ul>
 * Attempts to cast other Data types will produce an IllegalArgumentException by processDataInfo.
 * 
 * @author David Kofke
 *  
 */

/*
 * History Created on August 20, 2005 by kofke
 */
public class CastToInteger extends DataProcessor {

    /**
     * Sole constructor.
     */
    public CastToInteger() {
    }

    /**
     * Prepares processor to handle Data. Uses given DataInfo to determine the
     * type of Data to expect in subsequent calls to processData.
     * 
     * @throws IllegalArgumentException
     *             if input Data type is not one of those described in the
     *             general class comments
     */
    protected DataInfo processDataInfo(DataInfo inputDataInfo) {
        dataInteger = new DataInteger();
        if (inputDataInfo instanceof DataInfoDouble) {
            inputType = 0;
        } else if (inputDataInfo instanceof DataInfoDoubleArray) {
            inputType = 1;
        } else if (inputDataInfo instanceof DataInfoInteger) {
            inputType = 2;
            dataInteger = null;
        } else {
            throw new IllegalArgumentException("Cannot cast to int from "+ inputDataInfo.getClass());
        }
        DataInfo outputDataInfo = new DataInfoInteger(inputDataInfo.getLabel(), inputDataInfo.getDimension());
        outputDataInfo.addTags(inputDataInfo.getTags());
        return outputDataInfo;
    }
    
    /**
     * Extracts an int from the input data and returns it encapsulated in a
     * DataInteger.
     * 
     * @param data
     *            a Data instance of the type indicated by the DataInfo at
     *            construction
     * @return a DataInteger holding the value cast from the given Data; the same
     *         instance is returned with every invocation.
     * 
     * @throws ClassCastException
     *             if the given Data is not of the same type as indicated by the
     *             DataInfo given at construction
     */
    protected Data processData(Data data) {
        switch (inputType) {
        case 0:
            //we don't add ourselves
            dataInteger.x = (int)((DataDouble) data).x;
            return dataInteger;
        case 1:
            //we don't add ourselves
            dataInteger.x = (int)((DataDoubleArray) data).getData()[0];
            return dataInteger;
        case 2:
            return data;
        default:
            throw new Error("Assertion error.  Input type out of range: "+inputType);
        }
    }
    
    /**
     * Returns null.
     */
    public DataProcessor getDataCaster(DataInfo info) {
        return null;
    }

    private DataInteger dataInteger;
    private int inputType;
}
