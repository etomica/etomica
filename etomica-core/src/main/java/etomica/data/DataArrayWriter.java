/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.data;

import etomica.data.DataLogger.DataWriter;
import etomica.data.types.CastGroupToDoubleArray;
import etomica.data.types.CastToDoubleArray;
import etomica.data.types.DataDoubleArray;
import etomica.data.types.DataDoubleArray.DataInfoDoubleArray;
import etomica.data.types.DataGroup.DataInfoGroup;

import java.io.FileWriter;
import java.io.IOException;

/**
 * A DataWriter that writes out data as a table with optional column headings
 */
public class DataArrayWriter implements DataWriter, java.io.Serializable {

    private static final long serialVersionUID = 1L;
    private FileWriter fileWriter;
    private boolean firstWrite;
    private IDataInfo dataInfo;
    private boolean includeHeader;
    
    public DataArrayWriter() {
        super();
        setIncludeHeader(true);
        reset();
    }
    
    /**
     * Returns whether the writer includes column headers or not
     */
    public boolean getIncludeHeader() {
        return includeHeader;
    }

    /**
     * Directs the writer to include column headers or not
     */
    public void setIncludeHeader(boolean newIncludeHeader) {
        includeHeader = newIncludeHeader;
    }

    public void putDataInfo(IDataInfo newDataInfo) {
        dataInfo = newDataInfo;
    }

    public DataPipe getDataCaster(IDataInfo newDataInfo) {
        if (newDataInfo instanceof DataInfoDoubleArray) {
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
            // turn group into a multi-dimensional array, which we'll cast to a table
            // with a CastToTable
            return new CastGroupToDoubleArray();
        }
        return new CastToDoubleArray();
    }

    public void putData(IData data) {
        DataDoubleArray dataArray = (DataDoubleArray)data;
        try {
            int dim = dataArray.getArrayDimension();
            if (dim == 2) {
                int nColumns = dataArray.getArrayShape(1);
                if (nColumns < 1) {
                    return;
                }
                if (firstWrite && includeHeader) {
                    // if this is the first write to a file, start with the column headers
                    fileWriter.write(dataInfo.getLabel() + " " + dataInfo.getDimension() + "\n");
                }
                firstWrite = false;
                int[] indices = new int[2];
                for (int i = 0; i < dataArray.getArrayShape(0); i++) {
                    indices[0] = i;
                    indices[1] = 0;
                    fileWriter.write(Double.toString(dataArray.getValue(indices)));
                    for (int j = 1; j < nColumns; j++) {
                        indices[1] = j;
                        fileWriter.write("  " + Double.toString(dataArray.getValue(indices)));
                    }
                    fileWriter.write("\n");
                }
            } else {
                if (firstWrite && includeHeader) {
                    // if this is the first write to a file, start with the column headers
                    fileWriter.write(dataInfo.getLabel() + " " + dataInfo.getDimension() + "\n");
                }
                firstWrite = false;
                int[] indices = new int[2];
                for (int i = 0; i < dataArray.getArrayShape(0); i++) {
                    fileWriter.write("  " + Double.toString(dataArray.getValue(i)));
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
}
