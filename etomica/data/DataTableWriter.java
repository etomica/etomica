package etomica.data;

import java.io.FileWriter;
import java.io.IOException;

import etomica.data.DataLogger.DataWriter;
import etomica.data.types.CastGroupOfTablesToDataTable;
import etomica.data.types.CastGroupToDoubleArray;
import etomica.data.types.CastToTable;
import etomica.data.types.DataDoubleArray;
import etomica.data.types.DataTable;
import etomica.data.types.DataGroup.DataInfoGroup;
import etomica.data.types.DataTable.DataInfoTable;

/**
 * A DataWriter that writes out data as a table with column headings
 */
public class DataTableWriter implements DataWriter, java.io.Serializable {

    public DataTableWriter() {
        super();
        firstWrite = true;
    }
    
    public void putDataInfo(DataInfo newDataInfo) {
        dataInfo = newDataInfo;
    }

    public DataProcessor getDataCaster(DataInfo newDataInfo) {
        if (newDataInfo instanceof DataInfoTable) {
            return null;
        }
        else if (newDataInfo instanceof DataInfoGroup) {
            DataInfo dataInfo0 = ((DataInfoGroup)newDataInfo).getSubDataInfo(0);
            for (int i = 1; i<((DataInfoGroup)newDataInfo).getNDataInfo(); i++) {
                DataInfo subDataInfo = ((DataInfoGroup)newDataInfo).getSubDataInfo(0);
                if (subDataInfo.getClass() != dataInfo0.getClass()){
//                        ((DataInfoGroup)newDataInfo).getSubDataInfo(0))) {
                    throw new IllegalArgumentException("DataSinkTable can only handle homogeneous groups");
                }
            }
            if(((DataInfoGroup)newDataInfo).getSubDataInfo(0) instanceof DataInfoTable) {
                return new CastGroupOfTablesToDataTable();
            }
            return new CastGroupToDoubleArray();
        }
        return new CastToTable();
    }
    
    public void putData(Data data) {
        DataTable table = (DataTable)data;
        try {
            int nColumns = table.getNData();
            if (nColumns < 1) {
                return;
            }
            if (firstWrite) {
                // if this is the first write to a file, start with the column headers
                fileWriter.write(dataInfo.getLabel()+" "+dataInfo.getDimension()+"\n");
                fileWriter.write(((DataInfoTable)dataInfo).getSubDataInfo(0).getLabel());
                for (int i=1; i<nColumns; i++) {
                    fileWriter.write("  "+((DataInfoTable)dataInfo).getSubDataInfo(i).getLabel());
                }
                fileWriter.write("\n");
                firstWrite = false;
            }
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

    private FileWriter fileWriter;
    private boolean firstWrite;
    private DataInfo dataInfo;
}
