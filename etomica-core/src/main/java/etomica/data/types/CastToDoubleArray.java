/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.data.types;

import etomica.data.DataProcessor;
import etomica.data.IData;
import etomica.data.IDataInfo;
import etomica.data.types.DataDouble.DataInfoDouble;
import etomica.data.types.DataDoubleArray.DataInfoDoubleArray;
import etomica.data.types.DataTensor.DataInfoTensor;
import etomica.data.types.DataVector.DataInfoVector;

/**
 * A DataProcessor that converts a Data instance into a DataDoubleArray. Copies
 * an element of the input data to the DataDouble's encapsulated value and
 * returns the DataDouble instance.
 * <p>
 * Casting for various types of Data is performed as follows:
 * <ul>
 * <li><u>DataDoubleArray</u>. Does nothing, and returns input data directly.
 * 
 * <li><u>DataDouble</u>. Casts to a 1-element array.
 * 
 * <li><u>DataInteger</u>. Casts to a 1-element array.
 * 
 * <li><u>DataVector</u>. Casts to a one-dimensional array of length equal to
 * the number of vector elements.
 * 
 * <li><u>DataTensor</u>. Casts to a two-dimension square array.
 * 
 * <li><u>DataFunction</u>. Casts only the dependent data, to an array of the
 * same shape. Independent data are discarded.
 * 
 * <li><u>DataGroup</u>. Uses DataExtractor to locate a DataDoubleArray in DataGroup.  Other 
 * Data in group are discarded. 
 * 
 * </ul>
 * Attempts to process a different type will result in an
 * IllegalArgumentException when encountered in the processDataInfo method.
 * <p>
 * @author David Kofke and Andrew Schultz
 *  
 */
public class CastToDoubleArray extends DataProcessor {

    private DataDoubleArray outputData;
    private int inputType;

    /**
     * Sole constructor.
     */
    public CastToDoubleArray() {
    }

    /**
     * Prepares processor to perform cast. Given DataInfo is examined to see
     * what data type will be given to processor.
     *
     * @throws IllegalArgumentException
     *             if DataInfo is not one of the acceptable types, as described
     *             in general comments for this class
     */
    protected IDataInfo processDataInfo(IDataInfo inputDataInfo) {
        int[] arrayShape;
        if (inputDataInfo instanceof DataInfoDoubleArray) {
            inputType = 0;
            return inputDataInfo;
        } else if (inputDataInfo instanceof DataInfoDouble) {
            inputType = 1;
            arrayShape = new int[]{1};
        } else if (inputDataInfo instanceof DataInfoVector) {
            inputType = 2;
            arrayShape = new int[]{((DataInfoVector)inputDataInfo).getSpace().D()};
        } else if (inputDataInfo instanceof DataInfoTensor) {
            inputType = 3;
            int D = ((DataInfoTensor)inputDataInfo).getSpace().D();
            arrayShape = new int[]{D,D};
        } else {
            throw new IllegalArgumentException(
                    "Cannot cast to DataDoubleArray from " + inputDataInfo.getClass());
        }
        outputData = new DataDoubleArray(arrayShape);
        DataInfoDoubleArray outputDataInfo = new DataInfoDoubleArray(inputDataInfo.getLabel(), inputDataInfo.getDimension(), arrayShape);
        outputDataInfo.addTags(inputDataInfo.getTags());
        return outputDataInfo;
    }

    /**
     * Copies input Data to a DataDoubleArray and returns it (the DataDataDoubleArray).
     *
     * @throws ClassCastException
     *             if input Data is not of the type indicated by the most recent
     *             call to processDataInfo
     *
     */
    protected IData processData(IData data) {
        switch (inputType) {
        case 0:
            return data;
        case 1:
            outputData.E(((DataDouble) data).x);
            break;
        case 2:
            ((DataVector) data).x.assignTo(outputData.getData());
            break;
        case 3:
            // both Tensor and DataDoubleArray sequence data by rows
            ((DataTensor) data).x.assignTo(outputData.getData());
            break;
        default:
            throw new Error("Assertion error.  Input type out of range: "
                    + inputType);
        }
        //we don't add ourselves
        return outputData;
    }
}
