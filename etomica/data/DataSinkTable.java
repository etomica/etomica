package etomica.data;

import java.io.Serializable;
import java.util.HashMap;

import etomica.data.types.CastArrayToDoubleArray;
import etomica.data.types.CastGroupOfTablesToDataTable;
import etomica.data.types.CastGroupToDoubleArray;
import etomica.data.types.CastToTable;
import etomica.data.types.DataArray;
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
public class DataSinkTable implements DataSink, Serializable {

    public DataSinkTable() {
    }
    
    /* (non-Javadoc)
     * @see etomica.DataSink#getDataCaster(etomica.DataInfo)
     */
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
        } else if(dataInfo.getDataClass() == DataArray.class) {
            return new CastArrayToDoubleArray();
        }
        return new CastToTable();
    }
    /* (non-Javadoc)
     * @see etomica.DataSink#putData(etomica.Data)
     */
    public void putData(Data data) {
        Integer indexObj = (Integer)columnIndexHash.get(data);
        if (indexObj == null) {
            addNewData((DataTable)data);
            indexObj = new Integer(casterIndex++);
            if ((casterIndex+31)/32 > casterChangedBits.length) {
                casterChangedBits = Arrays.resizeArray(casterChangedBits,casterChangedBits.length+1);
            }
            casterChangedBits[casterChangedBits.length-1] |= 1<<((casterIndex-1)%32);
            columnIndexHash.put(data,indexObj);
        }
        int index = indexObj.intValue();
        casterChangedBits[index >> 5] &= (0xFFFFFFFF ^ (1<<(index & 0xFFFFFFFF)));
        if (updatingOnAnyChange) {
            fireDataChangedEvent();
        }
        else {
            int l = casterChangedBits.length;
            for (int i=0; i<l; i++) {
                if (casterChangedBits[i] != 0) {
                    return;
                }
            }
            fireDataChangedEvent();
            for (int i=0; i<l-1; i++) {
                casterChangedBits[i] = 0xFFFFFFFF;
            }
            if (casterIndex%32 == 0) {
                casterChangedBits[l-1] = 0xFFFFFFFF;
            }
            else {
                casterChangedBits[l-1] = (1<<(casterIndex%32)) - 1;
            }
        }
    }

    protected void addNewData(DataTable data) {
        int nColumns = data.getNColumns();
        if (dataColumns.length == 0) {
            rowHeaders = new String[data.getNRows()];
            for (int i=0; i<rowHeaders.length; i++) {
                rowHeaders[i] = data.getRowHeaders(i);
            }
        }
        int oldSize = dataColumns.length;
        dataColumns = (DataTable.Column[])Arrays.resizeArray(dataColumns,oldSize+nColumns);
        for (int i=0; i<nColumns; i++) {
            dataColumns[oldSize+i] = data.getColumn(i);
        }
        fireColumnCountChangedEvent();
        updateRowCount();
    }
    
    /* (non-Javadoc)
     * @see etomica.DataSink#putDataInfo(etomica.DataInfo)
     */
    public void putDataInfo(DataInfo dataInfo) {
    }
    
    public void reset() {
        columnIndexHash.clear();
        dataColumns = new DataTable.Column[0];
        updateRowCount();
        casterIndex = 0;
        casterChangedBits = new int[0];
        fireColumnCountChangedEvent();
        rowHeaders = null;
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

    /**
     * @return flag defining protocol for firing table events that
     * notify of column data changes.
     */
    public boolean isUpdatingOnAnyChange() {
        return updatingOnAnyChange;
    }

    /**
     * Describes protocol for firing table events that notify of column data
     * changes. If true, event is fired when any column notifies that it has
     * changed; if false, event is fired only when all columns have
     * notified that they have changed (since last time event was fired).
     */
    public void setUpdatingOnAnyChange(boolean updatingWithAnyChange) {
        this.updatingOnAnyChange = updatingWithAnyChange;
    }

    protected void fireDataChangedEvent() {
        for (int i = 0; i < listeners.length; i++) {
            listeners[i].tableDataChanged(this);
        }
    }

    protected void fireColumnCountChangedEvent() {
        for (int i = 0; i < listeners.length; i++) {
            listeners[i].tableColumnCountChanged(this);
        }
    }
    protected void fireRowCountChangedEvent() {
        for (int i = 0; i < listeners.length; i++) {
            listeners[i].tableRowCountChanged(this);
        }
    }
    
    public void addTableListener(DataTableListener newListener) {
        listeners = (DataTableListener[]) Arrays.addObject(listeners,
                newListener);
    }

    public void removeTableListener(DataTableListener oldListener) {
        listeners = (DataTableListener[]) Arrays.removeObject(listeners,
                oldListener);
    }

    private void updateRowCount() {
        int oldRowCount = rowCount;
        rowCount = 0;
        for (int i = 0; i < dataColumns.length; i++) {
            int n = dataColumns[i].getData().length;
            if (n > rowCount) {
                rowCount = n;
            }
        }
        if (oldRowCount != rowCount) {
            fireRowCountChangedEvent();
        }
    }
    
    private DataTableListener[] listeners = new DataTableListener[0];
    private boolean updatingOnAnyChange;
    private int rowCount;
    private DataTable.Column[] dataColumns = new DataTable.Column[0];
    private int casterIndex;
    private int[] casterChangedBits = new int[0];
    private HashMap columnIndexHash = new HashMap();
    private String[] rowHeaders;
}
