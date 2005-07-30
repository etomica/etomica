package etomica.log;

import java.io.FileWriter;
import java.io.IOException;

import etomica.Data;
import etomica.DataInfo;
import etomica.DataSink;
import etomica.data.DataProcessor;
import etomica.data.types.CastToDouble;
import etomica.data.types.CastToDoubleArray;
import etomica.data.types.DataDouble;
import etomica.data.types.DataDoubleArray;
import etomica.data.types.DataGroup;
import etomica.data.types.DataInteger;


/**
 * Logs data to a file as it comes in.
 */
public class DataLogger implements DataSink, java.io.Serializable {

	private final String fileName;
	
	public DataLogger(String aFileName) {
		fileName = aFileName;
        try { 
            FileWriter fileWriter = new FileWriter(fileName,false);
            fileWriter.close();
        }catch(IOException e) {
            System.err.println("Cannot open "+fileName+", caught IOException: " + e.getMessage());
        }
	}

    public void putDataInfo(DataInfo info) {}
    
    public DataProcessor getDataCaster(DataInfo info) {
        if (info.getDataClass() == DataDouble.class && info.getDataClass() == DataDoubleArray.class) {
            return null;
        }
        if (info.getDataClass() == DataInteger.class) {
            return new CastToDouble();
        }
        return new CastToDoubleArray();
    }
    
	public void putData(Data data) {
        try { 
            FileWriter fileWriter = new FileWriter(fileName,true);
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
    
}
