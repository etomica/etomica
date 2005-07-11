package etomica.log;

import java.io.FileWriter;
import java.io.IOException;

import etomica.Data;
import etomica.DataInfo;
import etomica.DataSink;
import etomica.data.DataGroup;
import etomica.data.types.DataDouble;
import etomica.data.types.DataDoubleArray;


/**
 * Logs data to a file as it comes in.
 */
public class DataLogger implements DataSink, java.io.Serializable {

	private final String fileName;
    private DataInfo dataInfo;
	
	public DataLogger(DataInfo info, String aFileName) {
        dataInfo = info;
		fileName = aFileName;
        try { 
            FileWriter fileWriter = new FileWriter(fileName,false);
            fileWriter.close();
        }catch(IOException e) {
            System.err.println("Cannot open "+fileName+", caught IOException: " + e.getMessage());
        }
	}

	public void putData(Data data) {
        dataInfo = data.getDataInfo();
        try { 
            FileWriter fileWriter = new FileWriter(fileName,true);
            if (data instanceof DataGroup) {
                for(int i=0; i<((DataGroup)data).getNData(); i++) {
                    putData(fileWriter,((DataGroup)data).getData(i));
                }
            }
            if (data instanceof DataDoubleArray) {
                double[] values = ((DataDoubleArray)data).getData();
                for(int j=0; j<values.length; j++) {fileWriter.write(values[j]+" ");}
            }
            else if (data instanceof DataDouble) {
                fileWriter.write(Double.toString(((DataDouble)data).x));
            }
    		fileWriter.write("\n");
            fileWriter.close();
        } catch(IOException e) {
            System.err.println("Cannot writing to "+fileName+", caught IOException: " + e.getMessage());
        }
	}
    
    /**
     * @exception IOException
     * @param fileWriter
     * @param data
     */
    private void putData(FileWriter fileWriter, Data data) throws IOException {
        if (data instanceof DataGroup) {
            for(int i=0; i<((DataGroup)data).getNData(); i++) {
                putData(fileWriter,((DataGroup)data).getData(i));
            }
        }
        if (data instanceof DataDoubleArray) {
            double[] values = ((DataDoubleArray)data).getData();
            for(int j=0; j<values.length; j++) {fileWriter.write(values[j]+" ");}
        }
        else if (data instanceof DataDouble) {
            fileWriter.write(Double.toString(((DataDouble)data).x));
        }
        fileWriter.write("\n");
        fileWriter.close();
    }        
    
}
