/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.data;

import etomica.data.types.*;
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
public class DataSinkTable extends DataSet {

    public DataSinkTable() {
        super();
        backwardDataMap = new int[0];
        forwardDataMap = new int[0];
    }

    public IDataSink makeDataSink(IDataInfo inputDataInfo) {
        DataSetSink sink = super.makeDataSink();
        if (inputDataInfo instanceof DataInfoTable) {
            return sink;
        }
        if (inputDataInfo instanceof DataInfoGroup) {
            for (int i = 1; i < ((DataInfoGroup) inputDataInfo).getNDataInfo(); i++) {
                if (((DataInfoGroup) inputDataInfo).getSubDataInfo(i).getClass() != ((DataInfoGroup) inputDataInfo).getSubDataInfo(0).getClass()) {
                    throw new IllegalArgumentException("DataSinkTable can only handle homoegeneous groups");
                }
            }
            if (((DataInfoGroup) inputDataInfo).getSubDataInfo(0) instanceof DataInfoTable) {
                CastGroupOfTablesToDataTable caster = new CastGroupOfTablesToDataTable();
                caster.setDataSink(sink);
                return caster;
            }
            CastGroupToDoubleArray caster0 = new CastGroupToDoubleArray();
            CastToTable caster = new CastToTable();
            caster0.setDataSink(caster);
            caster.setDataSink(sink);
            return caster0;
        }
        CastToTable caster = new CastToTable();
        caster.setDataSink(sink);
        return caster;
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
            for (int i = 0; i < listeners.length; i++) {
                if (listeners[i] instanceof DataTableListener) {
                    ((DataTableListener)listeners[i]).tableRowHeadersChanged(this);
                }
            }
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
            int n = ((DataInfoTable)psuedoSinks[i].getDataInfo()).getNRows();
            if (n > rowCount) {
                rowCount = n;
            }
        }
        if (oldRowCount != rowCount) {
            fireRowCountChangedEvent();
        }
    }
    
    private static final long serialVersionUID = 1L;
    private int rowCount;
}
