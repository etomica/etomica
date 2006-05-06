package etomica.data;

import java.io.Serializable;

import etomica.data.types.CastGroupOfTablesToDataTable;
import etomica.data.types.CastGroupToDoubleArray;
import etomica.data.types.CastToTable;
import etomica.data.types.DataDoubleArray;
import etomica.data.types.DataGroup;
import etomica.data.types.DataTable;
import etomica.data.types.DataGroup.DataInfoGroup;
import etomica.data.types.DataTable.DataInfoTable;

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
        backwardDataMap = new int[0];
        forwardDataMap = new int[0];
    }
    
    public String getRowHeader(int i) {
        return ((DataInfoTable)psuedoSinks[0].getDataInfo()).getRowHeader(i);
    }
    
    public double getValue(int row, int column) {
        int iTable = backwardDataMap[column];
        DataDoubleArray dataColumn = (DataDoubleArray)((DataTable)psuedoSinks[iTable].getData()).getData(column-forwardDataMap[iTable]);
        if(dataColumn.getLength() <= row) return Double.NaN;
        return dataColumn.getValue(row);
    }

    protected void dataInfoChanged(DataSetSink dataSetSink) {
        super.dataInfoChanged(dataSetSink);
        if (dataSetSink.index != -1) {
            updateRowCount();
        }
    }

    protected void dataChanged(DataSetSink dataSetSink) {
        boolean newData = (dataSetSink.index == -1);
        super.dataChanged(dataSetSink);
        if (newData) {
            updateRowCount();
        }
    }
    
    /**
     * Returns the length of the longest column.
     */
    public int getRowCount() {
        return rowCount;
    }
    
    protected void fireRowCountChangedEvent() {
        for (int i = 0; i < listeners.length; i++) {
            if (listeners[i] instanceof DataTableListener) {
                ((DataTableListener)listeners[i]).tableRowCountChanged(this);
            }
        }
    }
    
    private void updateRowCount() {
        int oldRowCount = rowCount;
        rowCount = 0;
        for (int i = 0; i < psuedoSinks.length; i++) {
            int n = ((DataTable)psuedoSinks[i].getData()).getNRows();
            if (n > rowCount) {
                rowCount = n;
            }
        }
        if (oldRowCount != rowCount) {
            fireRowCountChangedEvent();
        }
    }
    
    private int rowCount;
    
    protected static class DataCasterJudgeTable implements DataCasterJudge, Serializable {

        public DataProcessor getDataCaster(DataInfo dataInfo) {
            if (dataInfo instanceof DataInfoTable) {
                return null;
            } else if(dataInfo instanceof DataInfoGroup) {
                for (int i = 1; i<((DataInfoGroup)dataInfo).getNDataInfo(); i++) {
                    if (((DataInfoGroup)dataInfo).getSubDataInfo(i).getClass() != ((DataInfoGroup)dataInfo).getSubDataInfo(0).getClass()) {
                        throw new IllegalArgumentException("DataSinkTable can only handle homoegeneous groups");
                    }
                }
                if(((DataInfoGroup)dataInfo).getSubDataInfo(0) instanceof DataInfoTable) {
                    return new CastGroupOfTablesToDataTable();
                }
                return new CastGroupToDoubleArray();
            }
            return new CastToTable();
        }

    }
}
