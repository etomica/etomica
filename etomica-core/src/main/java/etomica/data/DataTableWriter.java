/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.data;

import java.io.FileWriter;
import java.io.IOException;

import etomica.data.DataLogger.DataWriter;
import etomica.data.types.CastGroupOfTablesToDataTable;
import etomica.data.types.CastGroupToDoubleArray;
import etomica.data.types.CastToTable;
import etomica.data.types.DataDoubleArray;
import etomica.data.types.DataGroup.DataInfoGroup;
import etomica.data.types.DataTable;
import etomica.data.types.DataTable.DataInfoTable;

/**
 * A DataWriter that writes out data as a table with optional column headings
 */
public class DataTableWriter implements DataWriter, java.io.Serializable {

    public DataTableWriter() {
        super();
        setIncludeHeader(true);
        reset();
    }

    /**
     * Directs the writer to include column headers or not
     */
    public void setIncludeHeader(boolean newIncludeHeader) {
        includeHeader = newIncludeHeader;
    }
    
    /**
     * Returns whether the writer includes column headers or not
     */
    public boolean getIncludeHeader() {
        return includeHeader;
    }
    
    public void putDataInfo(IDataInfo newDataInfo) {
        dataInfo = newDataInfo;
    }

    public DataPipe getDataCaster(IDataInfo newDataInfo) {
        if (newDataInfo instanceof DataInfoTable) {
            // we like tables
            return null;
        }
        else if (newDataInfo instanceof DataInfoGroup) {
            if (((DataInfoGroup)newDataInfo).getNDataInfo() == 0) {
                //it's empty, turn it into an empty array
                return new CastGroupToDoubleArray();
            }
            IDataInfo dataInfo0 = ((DataInfoGroup)newDataInfo).getSubDataInfo(0);
            for (int i = 1; i<((DataInfoGroup)newDataInfo).getNDataInfo(); i++) {
                IDataInfo subDataInfo = ((DataInfoGroup)newDataInfo).getSubDataInfo(0);
                if (subDataInfo.getClass() != dataInfo0.getClass()){
                    throw new IllegalArgumentException("DataSinkTable can only handle homogeneous groups");
                }
            }
            if (dataInfo0 instanceof DataInfoTable) {
                // group of tables, just flatten it
                return new CastGroupOfTablesToDataTable();
            }
            // turn group into a multi-dimensional array, which we'll cast to a table
            // with a CastToTable
            return new CastGroupToDoubleArray();
        }
        return new CastToTable();
    }
    
    public void putData(IData data) {
        DataTable table = (DataTable)data;
        try {
            int nColumns = table.getNData();
            if (nColumns < 1) {
                return;
            }
            if (firstWrite && includeHeader) {
                // if this is the first write to a file, start with the column headers
                fileWriter.write(dataInfo.getLabel()+" "+dataInfo.getDimension()+"\n");
                fileWriter.write(((DataInfoTable)dataInfo).getSubDataInfo(0).getLabel());
                for (int i=1; i<nColumns; i++) {
                    fileWriter.write("  "+((DataInfoTable)dataInfo).getSubDataInfo(i).getLabel());
                }
                fileWriter.write("\n");
            }
            firstWrite = false;
            for (int i=0; i<table.getNRows(); i++) {
                fileWriter.write(Double.toString(((DataDoubleArray)table.getData(0)).getData()[i]));
                for (int j=1; j<nColumns; j++) {
                    fileWriter.write("  "+Double.toString(((DataDoubleArray)table.getData(j)).getData()[i]));
                }
                fileWriter.write("\n");
            }
        }
        catch (IOException e) {
            throw new RuntimeException(e);
        }
        
    }
    
    public void setFileWriter(FileWriter newFileWriter) {
        fileWriter = newFileWriter;
    }
    
    public void reset() {
        firstWrite = true;
    }

    private static final long serialVersionUID = 2L;
    // DataLogger will give the fileWriter back to us when it actually writes
    private transient FileWriter fileWriter;
    private boolean firstWrite;
    private IDataInfo dataInfo;
    private boolean includeHeader;
}
