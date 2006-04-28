package etomica.data.types;


import etomica.data.Data;
import etomica.data.DataFactory;
import etomica.data.DataInfo;
import etomica.data.DataProcessor;
import etomica.data.types.DataGroup.DataInfoGroup;
import etomica.data.types.DataTable.DataInfoTable;

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
     * @throws IllegalArgumentException
     *             if DataInfo does not indicate a DataGroup, or if DataInfo
     *             indicates that expected DataGroup will not be homogeneous
     */
    protected DataInfo processDataInfo(DataInfo inputDataInfo) {
        if (inputDataInfo.getDataClass() != DataGroup.class) {
            throw new IllegalArgumentException("can only cast from DataGroup");
        }
        DataFactory[] elementFactory = ((DataGroup.Factory)inputDataInfo.getDataFactory()).getDataFactory();
        int nColumns = 0;
        int nRows = -1;
        for (int i = 0; i<((DataInfoGroup)inputDataInfo).getNDataInfo(); i++) {
            DataInfo elementDataInfo = ((DataInfoGroup)inputDataInfo).getSubDataInfo(i);
            if (elementDataInfo.getDataClass() != DataTable.class) {
                throw new IllegalArgumentException("this class can only cast groups of DataTables");
            }
            nColumns += ((DataInfoTable)elementDataInfo).getNDataInfo();
            if (nRows > -1 && ((DataInfoTable)elementDataInfo).getNRows() != nRows) {
                throw new IllegalArgumentException("all columns must have an equal number of rows");
            }
            nRows = ((DataTable.Factory)elementFactory[i]).getNRows();
        }
        
        outputData = new DataTable(nColumns, nRows);
        
        DataFactory outputFactory = DataTable.getFactory(nRows, nColumns);
        return new DataInfo(inputDataInfo.getLabel(), inputDataInfo.getDimension(), outputFactory);
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
        return outputData;
    }
    
    /**
     * Returns null.
     */
    public DataProcessor getDataCaster(DataInfo info) {
        if (info.getDataClass() != DataGroup.class) {
            throw new IllegalArgumentException("can only cast from DataGroup");
        }
        return null;
    }

    private DataTable outputData;
}
