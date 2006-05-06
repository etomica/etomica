package etomica.data.types;


import etomica.data.Data;
import etomica.data.DataInfo;
import etomica.data.DataProcessor;
import etomica.data.types.DataDoubleArray.DataInfoDoubleArray;
import etomica.data.types.DataGroup.DataInfoGroup;
import etomica.data.types.DataTable.DataInfoTable;
import etomica.util.Arrays;

/**
 * A DataProcessor that converts a homogeneous DataGroup into a multidimensional DataDoubleArray. 
 * All elements in group must be of the same data type.  If the group contains only one element,
 * the effect is the same as applying CastToDoubleArray to the element.  Otherwise, if <tt>d</tt>
 * is the dimension of each grouped element, the resulting DataDoubleArray will be of 
 * dimension <tt>d+1</tt>.  The added dimension corresponds to the different elements of the group.
 * <p>
 * For example:
 * <ul>
 * <li>if the DataGroup contains 8 DataDouble instances (thus d = 0), the cast gives a
 * one-dimensional DataDoubleArray of length 8.
 * <li>if the DataGroup contains 4 Vector3D instances (thus d = 1), the cast gives a two-dimensional DataDoubleArray
 * of dimensions (4, 3).
 * <li>if the DataGroup contains 5 two-dimensional (d = 2) DataDoubleArray of shape (7,80), the cast gives a
 * three-dimensional DataDoubleArray of shape (5, 7, 80).
 * </ul>
 * etc. 
 * 
 * @author David Kofke
 *  
 */

/*
 * History Created on Jul 21, 2005 by kofke
 */
public class CastGroupOfTablesToDataTable extends DataProcessor {

    /**
     * Sole constructor.
     */
    public CastGroupOfTablesToDataTable() {
    }

    /**
     * Prepares processor to handle Data. Uses given DataInfo to determine the
     * type of Data to expect in subsequent calls to processData.
     * 
     * @throws ClassCastException
     *             if DataInfo does not indicate a DataGroup or its sub-groups 
     *             are not DataTables
     *             
     * @throws IllegalArgumentException
     *             if DataInfo indicates that the DataTables have different
     *             numbers of rows
     */
    protected DataInfo processDataInfo(DataInfo inputDataInfo) {
        if (!(inputDataInfo instanceof DataInfoGroup)) {
            throw new IllegalArgumentException("can only cast from DataGroup");
        }
        DataInfoDoubleArray[] columnDataInfo = new DataInfoDoubleArray[0];
        int nColumns = 0;
        int nRows = -1;
        for (int i = 0; i<((DataInfoGroup)inputDataInfo).getNDataInfo(); i++) {
            DataInfoTable elementDataInfo = (DataInfoTable)((DataInfoGroup)inputDataInfo).getSubDataInfo(i);
            columnDataInfo = (DataInfoDoubleArray[])Arrays.resizeArray(columnDataInfo, nColumns+elementDataInfo.getNDataInfo());
            for (int j=nColumns; i<columnDataInfo.length; i++) {
                columnDataInfo[j] = (DataInfoDoubleArray)elementDataInfo.getSubDataInfo(j-nColumns);
            }
            nColumns = columnDataInfo.length;
            if (nRows > -1 && elementDataInfo.getNRows() != nRows) {
                throw new IllegalArgumentException("all columns must have an equal number of rows");
            }
            nRows = elementDataInfo.getNRows();
        }
        
        outputData = null;
        
        return new DataInfoTable(inputDataInfo.getLabel(), columnDataInfo, nRows, null);
    }
    
    /**
     * Converts data in given group to a DataDoubleArray as described in the
     * general comments for this class.
     * 
     * @throws ClassCastException
     *             if the given Data is not a DataGroup with Data elements of
     *             the type indicated by the most recent call to
     *             processDataInfo.
     */
    protected Data processData(Data data) {
        if (outputData == null) {
            DataDoubleArray[] columns = new DataDoubleArray[outputDataInfo.getNDataInfo()];
            int i=0;
            for (int j=0; j<((DataGroup)data).getNData(); j++) {
                for (int k=0; k<((DataTable)((DataGroup)data).getData(j)).getNData(); k++) {
                    columns[i] = (DataDoubleArray)((DataTable)((DataGroup)data).getData(j)).getData(i);
                }
            }
            outputData = new DataTable(columns);
        }
        return outputData;
    }
    
    /**
     * Returns null.
     */
    public DataProcessor getDataCaster(DataInfo info) {
        if (!(info instanceof DataInfoGroup)) {
            throw new IllegalArgumentException("can only cast from DataGroup");
        }
        return null;
    }

    private DataTable outputData;
    private DataInfoTable outputDataInfo;
}
