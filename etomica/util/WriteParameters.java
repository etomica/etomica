package etomica.util;

import java.io.FileWriter;
import java.io.IOException;
import java.lang.reflect.Field;

/**
 * This class writes input parameters from a parameter wrapper Object to a file.
 */
public class WriteParameters implements java.io.Serializable {

    public WriteParameters() {}

    public WriteParameters(String inputFileName) {
        this();
        setInputFileName(inputFileName);
    }
    
    public WriteParameters(String inputFileName, ParameterBase parameterWrapper) {
        this(inputFileName);
        setParameterWrapper(parameterWrapper);
    }
    
    /**
     * Returns the fileName parameters are read from.
     */
    public String getInputFileName() {
        return fileName;
    }

    /**
     * Sets the fileName to read the parameters from.
     */
    public void setInputFileName(String newFileName) {
        fileName = newFileName;
    }

    /**
     * Returns the parameter wrapper.
     */
    public ParameterBase getParameterWrapper() {
        return wrapper;
    }

    /**
     * Sets the parameterWrapper
     */
    public void setParameterWrapper(ParameterBase newParameterWrapper) {
        wrapper = newParameterWrapper;
    }

    /**
     * Reads each line from inputFileName and attempts to match it with a field
     * from parameterWrapper and sets the field to the value from the file.
     * This routine handles boolean, int, long, double, String and arrays of 
     * boolean, int long and double.
     */
    public void writeParameters() 
              throws IOException {
        if (fileName == null) {
            throw new IllegalStateException("you must set the file name to write");
        }
        if (wrapper == null) {
            throw new IllegalStateException("you must set the parameter wrapper object to write");
        }
        
        Field[] fields = wrapper.getClass().getFields();
        FileWriter fileWriter;
        try {
            fileWriter = new FileWriter(fileName);
        }
        catch(IOException e) {
            throw new RuntimeException("Cannot open "+fileName+", caught IOException: " + e.getMessage(),e);
        }

        fileWriter.write("# Class "+wrapper.getClass().getName()+"\n");
        // write fields that were not in the file before
        for (int j=0; j<fields.length; j++) {
            fileWriter.write(fields[j].getName()+" "+wrapper.getValueString(fields[j])+"\n");
        }
        fileWriter.close();
    }
    
    private static final long serialVersionUID = 1L;
    protected ParameterBase wrapper;
    protected String fileName;
}
