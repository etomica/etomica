package etomica.log;

import java.io.FileWriter;
import java.io.IOException;

import etomica.DataSink;
import etomica.units.Dimension;


/**
 * Logs data to a file as it comes in.
 */
public class DataLogger implements DataSink {

	private final String fileName;
    private Dimension dimension;
    private String label;
	
	public DataLogger(String aFileName) {
		fileName = aFileName;
        try { 
            FileWriter fileWriter = new FileWriter(fileName,false);
            fileWriter.close();
        }catch(IOException e) {
            System.err.println("Cannot open "+fileName+", caught IOException: " + e.getMessage());
        }
	}

	public void putData(double[] values) {
        try { 
            FileWriter fileWriter = new FileWriter(fileName,true);
    		for(int j=0; j<values.length; j++) {fileWriter.write(values[j]+" ");}
    		fileWriter.write("\n");
            fileWriter.close();
        }catch(IOException e) {
            System.err.println("Cannot writing to "+fileName+", caught IOException: " + e.getMessage());
        }
	}
    
    public void setDimension(Dimension dimension) {
        this.dimension = dimension;
    }

    public void setLabel(String label) {
        this.label = label;
    }
}
