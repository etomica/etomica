package etomica.log;

import java.io.FileWriter;
import java.io.IOException;

import etomica.DataSink;


/**
 * Logs data to a file as it comes in.
 */
public class DataLogger implements DataSink {

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

	public void add(double[] values) {
        try { 
            FileWriter fileWriter = new FileWriter(fileName,true);
    		for(int j=0; j<values.length; j++) {fileWriter.write(values[j]+" ");}
    		fileWriter.write("\n");
            fileWriter.close();
        }catch(IOException e) {
            System.err.println("Cannot writing to "+fileName+", caught IOException: " + e.getMessage());
        }
	}

}
