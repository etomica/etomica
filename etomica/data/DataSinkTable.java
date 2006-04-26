package etomica.data;

import java.io.Serializable;

import etomica.data.types.CastGroupOfTablesToDataTable;
import etomica.data.types.CastGroupToDoubleArray;
import etomica.data.types.CastToTable;
import etomica.data.types.DataGroup;
import etomica.data.types.DataTable;
import etomica.util.Arrays;

/**
 * Receives data from multiple streams and organizes it in the form of a table,
 * with each stream of data forming a column of the table. Fires event when data
 * in table is changed; can be configured to fire when any column of data
 * changes, or to perform action only when all columns are changed.
 * 
 * @author David Kofke
 *  
 */

/*
 * History Created on Apr 9, 2005 by kofke
 */
public class DataSinkTable extends DataSet {

    public DataSinkTable() {
        super(new DataCasterJudgeTable());
    }
    
    /* (non-Javadoc)
     * @see etomica.DataSink#getDataCaster(etomica.DataInfo)
     */

    protected void addNewData(Data data) {
        DataTable dataTable = (DataTable)data;
        int nColumns = dataTable.getNColumns();
        if (dataColumns.length == 0) {
            rowHeaders = new String[dataTable.getNRows()];
            for (int i=0; i<rowHeaders.length; i++) {
                rowHeaders[i] = dataTable.getRowHeaders(i);
            }
        }
        int oldSize = dataColumns.length;
        dataColumns = (DataTable.Column[])Arrays.resizeArray(dataColumns,oldSize+nColumns);
        for (int i=0; i<nColumns; i++) {
            dataColumns[oldSize+i] = dataTable.getColumn(i);
        }
        super.addNewData(data);
        updateRowCount(dataTable);
    }
    
    public void reset() {
        dataColumns = new DataTable.Column[0];
        super.reset();
        updateRowCount(null);
    }
    
    
    public double getValue(int row, int column) {
        if(dataColumns[column].getData().length <= row) return Double.NaN;
        return dataColumns[column].getData()[row];
    }

    public int getColumnCount() {
        return dataColumns.length;
    }

    /**
     * Returns the length of the longest column.
     */
    public int getRowCount() {
        return rowCount;
    }
    
    public String getRowHeader(int i) {
        return rowHeaders[i];
    }

    public DataTable.Column getColumn(int column) {
        return dataColumns[column];
    }
    
//    /**
//     * Descriptive string for each column.  Returns an
//     * array of length equal to the number of columns in
//     * the table.  Each entry is obtained from the getLabel
//     * method for the column. Never returns null.
//     */
//    public String[] getColumnHeadings() {
//        //need to refresh with each call because columnHeadings isn't 
//        //notified when setLabel of DataBin is called
//        for (int i = 0; i < columns.length; i++) {
//            columnHeadings[i] = columns[i].getDataInfo().getLabel();
//        }
//        return columnHeadings;
//    }

    protected void fireRowCountChangedEvent() {
        for (int i = 0; i < listeners.length; i++) {
            if (listeners[i] instanceof DataTableListener) {
                ((DataTableListener)listeners[i]).tableRowCountChanged(this);
            }
        }
    }
    
    private void updateRowCount(DataTable dataTable) {
        int oldRowCount = rowCount;
        rowCount = 0;
        for (int i = 0; i < dataColumns.length; i++) {
            int n = dataColumns[i].getData().length;
            if (n > rowCount) {
                rowCount = n;
            }
        }
        if (oldRowCount != rowCount) {
            if (dataTable == null) {
                rowHeaders = new String[0];
            }
            else {
                rowHeaders = new String[dataTable.getNRows()];
                for (int i=0; i<rowHeaders.length; i++) {
                    rowHeaders[i] = dataTable.getRowHeaders(i);
                }
            }
            fireRowCountChangedEvent();
        }
    }
    
    private int rowCount;
    private DataTable.Column[] dataColumns = new DataTable.Column[0];
    private String[] rowHeaders;
    
    protected static class DataCasterJudgeTable implements DataCasterJudge, Serializable {

        public DataProcessor getDataCaster(DataInfo dataInfo) {
            if (dataInfo.getDataClass() == DataTable.class) {
                return null;
            } else if(dataInfo.getDataClass() == DataGroup.class) {
                DataInfo[] info = ((DataGroup.Factory)dataInfo.getDataFactory()).getDataInfoArray(); 
                Class innerDataClass = info[0].getDataClass();
                for (int i = 1; i<info.length; i++) {
                    if (info[i].getDataClass() != innerDataClass) {
                        throw new IllegalArgumentException("DataSinkTable can only handle homogeneous groups");
                    }
                }
                if(innerDataClass == DataTable.class) {
                    return new CastGroupOfTablesToDataTable();
                }
                return new CastGroupToDoubleArray();
            }
            return new CastToTable();
        }

    }
}
