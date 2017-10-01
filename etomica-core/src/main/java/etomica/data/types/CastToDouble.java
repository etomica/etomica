/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.data.types;

import etomica.data.DataPipe;
import etomica.data.DataProcessor;
import etomica.data.IData;
import etomica.data.IDataInfo;
import etomica.data.types.DataDouble.DataInfoDouble;
import etomica.data.types.DataDoubleArray.DataInfoDoubleArray;

/**
 * A DataProcessor that converts a Data instance into a DataDouble. Copies an
 * element of the input data to the DataDouble's encapsulated value and returns
 * the DataDouble instance.
 * Casting for various types of Data is performed as follows:
 * <ul>
 * <li><u>DataDoubleArray</u>. Uses only first element of array.
 * 
 * <li><u>DataDouble</u>. Does nothing, and returns input data directly.
 * 
 * <li><u>DataInteger</u>. Performs a simple cast of the int to a double.
 * 
 * <li><u>DataGroup</u>. Uses DataExtractor to locate a DataDouble in DataGroup.  Other 
 * Data in group are discarded. 
 * </ul>
 * Attempts to cast other Data types will produce an IllegalArgumentException by processDataInfo.
 * <p>
 * @author David Kofke
 *  
 */
public class CastToDouble extends DataProcessor {

    /**
     * Sole constructor.
     */
    public CastToDouble() {
    }

    /**
     * Prepares processor to handle Data. Uses given DataInfo to determine the
     * type of Data to expect in subsequent calls to processData.
     * 
     * @throws IllegalArgumentException
     *             if input Data type is not one of those described in the
     *             general class comments
     */
    protected IDataInfo processDataInfo(IDataInfo inputDataInfo) {
        dataDouble = new DataDouble();
        if (inputDataInfo instanceof DataInfoDouble) {
            inputType = 0;
        } else if (inputDataInfo instanceof DataInfoDoubleArray) {
            inputType = 1;
        } else {
            throw new IllegalArgumentException("Cannot cast to double from "
                    + inputDataInfo.getClass());
        }
        
        DataInfoDouble outputDataInfo = new DataInfoDouble(inputDataInfo.getLabel(), inputDataInfo.getDimension());
        outputDataInfo.addTags(inputDataInfo.getTags());
        return outputDataInfo;
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
    protected IData processData(IData data) {
        switch (inputType) {
        case 0:
            return data;
        case 1:
            //we don't add ourselves
            dataDouble.x = ((DataDoubleArray) data).getData()[0];
            return dataDouble;
        default:
            throw new Error("Assertion error.  Input type out of range: "+inputType);
        }
    }
    
    /**
     * Returns null.
     */
    public DataPipe getDataCaster(IDataInfo info) {
        return null;
    }

    private static final long serialVersionUID = 1L;
    private DataDouble dataDouble;
    private int inputType;
}
