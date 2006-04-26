package etomica.data;

import java.io.FileWriter;
import java.io.IOException;

import etomica.data.DataLogger.DataWriter;
import etomica.data.types.CastGroupOfTablesToDataTable;
import etomica.data.types.CastGroupToDoubleArray;
import etomica.data.types.CastToTable;
import etomica.data.types.DataGroup;
import etomica.data.types.DataTable;

/**
 * A DataWriter that writes out data as a table with column headings
 */
public class DataTableWriter implements DataWriter, java.io.Serializable {

    public DataTableWriter() {
        super();
        firstWrite = true;
    }
    
    public void putDataInfo(DataInfo newDataInfo) {
        // no DataSink
    }

    public DataProcessor getDataCaster(DataInfo newDataInfo) {
        if (newDataInfo.getDataClass() == DataTable.class) {
            return null;
        }
        else if (newDataInfo.getDataClass() == DataGroup.class) {
            DataInfo[] info = ((DataGroup.Factory)newDataInfo.getDataFactory()).getDataInfoArray(); 
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
    
    public void putData(Data data) {
        DataTable table = (DataTable)data;
        try {
            int nColumns = table.getNColumns();
            if (nColumns < 1) {
                return;
            }
            if (firstWrite) {
                // if this is the first write to a file, start with the column headers
                fileWriter.write(table.getDataInfo().getLabel()+" "+table.getDataInfo().getDimension()+"\n");
                fileWriter.write(table.getColumn(0).getHeading());
                for (int i=1; i<nColumns; i++) {
                    fileWriter.write("  "+table.getColumn(i).getHeading());
                }
                fileWriter.write("\n");
                firstWrite = false;
            }
            for (int i=0; i<table.getNRows(); i++) {
                fileWriter.write(Double.toString(table.getColumn(0).getData()[i]));
                for (int j=1; j<nColumns; j++) {
                    fileWriter.write("  "+Double.toString(table.getColumn(j).getData()[i]));
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
}
