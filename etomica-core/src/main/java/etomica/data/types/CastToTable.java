/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.data.types;

import etomica.data.DataProcessor;
import etomica.data.DataSourceIndependent;
import etomica.data.IData;
import etomica.data.IDataInfo;
import etomica.data.types.DataDouble.DataInfoDouble;
import etomica.data.types.DataDoubleArray.DataInfoDoubleArray;
import etomica.data.types.DataFunction.DataInfoFunction;
import etomica.data.types.DataGroup.DataInfoGroup;
import etomica.data.types.DataTable.DataInfoTable;
import etomica.data.types.DataTensor.DataInfoTensor;
import etomica.data.types.DataVector.DataInfoVector;
import etomica.space.Tensor;

import java.io.Serializable;

/**
 * 
 * A DataProcessor that converts a Data instance into a DataTable. Casting for various
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
public class CastToTable extends DataProcessor implements Serializable {

    private DataTable outputData;
    private int inputType;
    private DataSourceIndependent xDataSource;

    /**
     * Sole constructor.
     */
    public CastToTable() {
    }

    /**
     * Prepares processor to perform cast. Given DataInfo is examined to see
     * what data type will be given to processor.
     *
     * @throws IllegalArgumentException
     *             if DataInfo is not one of the acceptable types, as described
     *             in general comments for this class
     */
    public IDataInfo processDataInfo(IDataInfo inputDataInfo) {
        if (inputDataInfo instanceof DataInfoGroup) {
            throw new IllegalArgumentException("Cannot cast to DataTable from "
                    + inputDataInfo.getClass());
        }
        int nRows, nColumns;
        String[] rowHeaders = null;
        IDataInfo[] columnInfo = new IDataInfo[]{inputDataInfo};
        if (inputDataInfo instanceof DataInfoFunction) {
            int[] arrayShape = ((DataInfoFunction)inputDataInfo).getArrayShape();
            if (arrayShape.length != 1) {
                throw new IllegalArgumentException("DataFunction must be 1-dimensional");
            }
            nColumns = 2;
            nRows = arrayShape[0];
            columnInfo = new IDataInfo[2];
            xDataSource = ((DataInfoFunction)inputDataInfo).getXDataSource();
            columnInfo[0] = xDataSource.getIndependentDataInfo(0);
            columnInfo[1] = inputDataInfo;
            inputType = 5;
        }
        else if (inputDataInfo instanceof DataInfoDoubleArray) {
            int[] arrayShape = ((DataInfoDoubleArray) inputDataInfo)
                    .getArrayShape();
            if (arrayShape.length > 2) {
                throw new IllegalArgumentException(
                        "Cannot cast to table a data set with dimension greater than 2: "
                                + arrayShape.length);
            }
            if (arrayShape.length == 1) {
                inputType = 0;
                nColumns = 1;
                nRows = arrayShape[0];
            } else {
                inputType = 1;
                nColumns = arrayShape[0];
                nRows = arrayShape[1];
                columnInfo = new IDataInfo[nColumns];
                for (int i=0; i<nColumns; i++) {
                    columnInfo[i] = inputDataInfo;
                }
            }
        } else if (inputDataInfo instanceof DataInfoDouble) {
            inputType = 2;
            nColumns = 1;
            nRows = 1;
        } else if (inputDataInfo instanceof DataInfoVector) {
            inputType = 3;
            nColumns = 1;
            nRows = ((DataInfoVector)inputDataInfo).getSpace().D();
        } else if (inputDataInfo instanceof DataInfoTensor) {
            inputType = 4;
            int D = ((DataInfoTensor)inputDataInfo).getSpace().D();
            nColumns = D;
            nRows = D;
            columnInfo = new IDataInfo[nColumns];
            for (int i=0; i<nColumns; i++) {
                columnInfo[i] = inputDataInfo;
            }
        } else {
            throw new IllegalArgumentException("Cannot cast to DataTable from "
                    + inputDataInfo.getClass());
        }
        outputData = new DataTable(nColumns, nRows);
        if (inputType == 0 || inputType == 5) {
            // we'll actually wrap the incoming DataDoubleArray(s) in a DataTable
            outputData = null;
        }
        return new DataInfoTable("Table", columnInfo, nRows, rowHeaders);

    }

    /**
     * Copies input Data to a DataTable and returns it (the DataTable).
     *
     * @throws ClassCastException
     *             if input Data is not of the type indicated by the most recent
     *             call to processDataInfo
     *
     */
    protected IData processData(IData data) {
        switch (inputType) {
        case 0: //DataDoubleArray
            if (outputData == null) {
                outputData = new DataTable(new DataDoubleArray[]{(DataDoubleArray)data});
            }
        case 1: //2D DataDoubleArray
            int nColumns = outputData.getNData();
            for (int i = 0; i < nColumns; i++) {
                ((DataDoubleArray) data).assignColumnTo(i,
                        ((DataDoubleArray)outputData.getData(i)).getData());
            }
            break;
        case 2: //DataDouble
            outputData.E(((DataDouble) data).x);
            break;
        case 3: //DataVector
            ((DataVector) data).x.assignTo(((DataDoubleArray)outputData.getData(0)).getData());
            break;
        case 4: //DataTensor
            Tensor x = ((DataTensor) data).x;
            for (int i = 0; i < x.D(); i++) {
                for (int j = 0; j < x.D(); j++) {
                    ((DataDoubleArray)outputData.getData(i)).getData()[j] = x.component(i, j);
                }
            }
            break;
        case 5: //DataFunction
            if (outputData == null) {
                outputData = new DataTable(new DataDoubleArray[]{xDataSource.getIndependentData(0),(DataDoubleArray)data});
            }
            break;
        }
        return outputData;
    }
}
