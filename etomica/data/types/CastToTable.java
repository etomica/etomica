package etomica.data.types;

import java.io.Serializable;

import etomica.data.Data;
import etomica.data.DataFactory;
import etomica.data.DataInfo;
import etomica.data.DataProcessor;
import etomica.space.Tensor;

/**
 * 
 * Processor that takes Data and converts it to a DataTable. Casting for various
 * types of Data is performed as follows:
 * <ul>
 * <li><u>DataDoubleArray</u>. A one-dimensional array will be put into a
 * single column, and a two-dimensional array will be arranged into a
 * corresponding set of multiple columns; attempts to cast higher-dimensional
 * arrays will results in an IllegalArgumentException.
 * 
 * <li><u>DataDouble</u>. Casts to a 1-column, 1-row table.
 * 
 * <li><u>DataInteger</u>. Casts to a 1-column, 1-row table.
 * 
 * <li><u>DataVector</u>. Places vector values in a single column.
 * 
 * <li><u>DataTensor</u>. Arranges elements into multiple columns, like a
 * matrix.
 * 
 * <li><u>DataFunction</u>. Handles only a 1-dimensional DataFunction. Casts
 * to a table of two columns, containing the independent and dependent data,
 * respectively.
 * </ul>
 * Attempts to process a different type will result in an
 * IllegalArgumentException when encountered in the processDataInfo method.
 * <p>
 * @author Andrew Schultz and David Kofke
 *  
 */

/*
 * History Created on Jul 28, 2005 by kofke
 */
public class CastToTable extends DataProcessor implements Serializable {

    /**
     * Prepares processor to perform cast. Given DataInfo is examined to see
     * what data type will be given to processor.
     * 
     * @throws IllegalArgumentException
     *             if DataInfo is not one of the acceptable types, as described
     *             in general comments for this class
     */
    public DataInfo processDataInfo(DataInfo inputDataInfo) {
        Class inputClass = inputDataInfo.getDataClass();
        DataFactory factory = inputDataInfo.getDataFactory();
        if (inputClass == DataDoubleArray.class) {
            inputType = 0;
            int[] arrayShape = ((DataDoubleArray.Factory) factory)
                    .getArrayShape();
            if (arrayShape.length > 2) {
                throw new IllegalArgumentException(
                        "Cannot cast to table a data set with dimension greater than 2: "
                                + arrayShape.length);
            }
            if (arrayShape.length == 1) {
                outputData = new DataTable(inputDataInfo.getLabel(),
                        new DataInfo[] { inputDataInfo }, arrayShape[0]);
            } else {
                outputData = new DataTable(inputDataInfo.getLabel(),
                        inputDataInfo.getDimension(), arrayShape[0],
                        arrayShape[1]);
            }
        } else if (inputClass == DataDouble.class) {
            inputType = 1;
            outputData = new DataTable(inputDataInfo.getLabel(),
                    new DataInfo[] { inputDataInfo }, 1);
        } else if (inputClass == DataInteger.class) {
            inputType = 2;
            outputData = new DataTable(inputDataInfo.getLabel(),
                    new DataInfo[] { inputDataInfo }, 1);
        } else if (inputClass == DataVector.class) {
            inputType = 3;
            int D = ((DataVector.Factory) factory).getSpace().D();
            outputData = new DataTable(inputDataInfo.getLabel(),
                    new DataInfo[] { inputDataInfo }, D);
        } else if (inputClass == DataTensor.class) {
            inputType = 4;
            int D = ((DataTensor.Factory) factory).getSpace().D();
            outputData = new DataTable(inputDataInfo.getLabel(), inputDataInfo
                    .getDimension(), D, D);
        } else if (inputClass == DataFunction.class) {
            int[] sizes = ((DataFunction.Factory) factory)
                    .getIndependentDataSizes();
            DataDoubleArray[] data = ((DataFunction.Factory) factory)
                    .getIndependentData();
            if (sizes.length == 1) {
                outputData = new DataTable(
                        inputDataInfo.getLabel(),
                        new DataInfo[] { data[0].getDataInfo(), inputDataInfo },
                        sizes[0]);
            } else {//if (sizes.length == 2) {
                throw new IllegalArgumentException(
                        "DataFunction must be 1-dimensional");
            }
            inputType = 5;
        } else {
            throw new IllegalArgumentException("Cannot cast to DataTable from "
                    + inputClass);
        }
        return outputData.getDataInfo();

    }//end of processDataInfo

    /**
     * Copies input Data to a DataTable and returns it (the DataTable).
     * 
     * @throws ClassCastException
     *             if input Data is not of the type indicated by the most recent
     *             call to processDataInfo
     *  
     */
    protected Data processData(Data data) {
        switch (inputType) {
        case 0: //DataDoubleArray
            int nColumns = outputData.myColumns.length;
            for (int i = 0; i < nColumns; i++) {
                ((DataDoubleArray) data).assignColumnTo(i,
                        outputData.myColumns[i].getData());
            }
            break;
        case 1: //DataDouble
            outputData.myColumns[0].getData()[0] = ((DataDouble) data).x;
            break;
        case 2: //DataInteger
            outputData.myColumns[0].getData()[0] = ((DataInteger) data).x;
            break;
        case 3: //DataVector
            ((DataVector) data).x.assignTo(outputData.myColumns[0].getData());
            break;
        case 4: //DataTensor
            Tensor x = ((DataTensor) data).x;
            for (int i = 0; i < x.length(); i++) {
                for (int j = 0; j < x.length(); j++) {
                    outputData.myColumns[i].getData()[j] = x.component(i, j);
                }
            }
            break;
        case 5: //DataFunction
            ((DataFunction) data).getXData(0).assignColumnTo(0,
                    outputData.myColumns[0].getData());
            ((DataFunction) data).getYData().assignColumnTo(0,
                    outputData.myColumns[1].getData());
            break;
        }
        return outputData;
    }

    /**
     * Returns null, indicating the this DataProcessor can handle (almost) any
     * Data type.
     */
    public DataProcessor getDataCaster(DataInfo dataInfo) {
        return null;
    }

    private DataTable outputData;
    private int inputType;
}