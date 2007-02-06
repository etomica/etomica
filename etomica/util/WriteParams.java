package etomica.util;

import java.io.FileWriter;
import java.io.IOException;
import java.lang.reflect.Field;

/**
 * This class writes input parameters from a parameter wrapper Object to a file.
 */
public class WriteParams implements java.io.Serializable {

    public WriteParams() {}

    public WriteParams(String inputFileName) {
        this();
        setInputFileName(inputFileName);
    }
    
    public WriteParams(String inputFileName, ParamBase parameterWrapper) {
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
    public ParamBase getParameterWrapper() {
        return wrapper;
    }

    /**
     * Sets the parameterWrapper
     */
    public void setParameterWrapper(ParamBase newParameterWrapper) {
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
            doWrite(fileWriter,fields[j]);
        }
        fileWriter.close();
    }
    
    /**
     * Actually writes field name and value to a file using fileWriter.
     * If the type of field is unrecognized, the wrapper object itself will
     * be asked to write the field.
     */
    protected void doWrite(FileWriter fileWriter, Field field) 
               throws IOException {
        Class type = field.getType();
        try {
            if (type == String.class) {
                fileWriter.write(field.getName()+" "+field.get(wrapper));
            }
            else if (type == int.class) {
                fileWriter.write(field.getName()+" "+field.getInt(wrapper));
            }
            else if (type == long.class) {
                fileWriter.write(field.getName()+" "+field.getLong(wrapper));
            }
            else if (type == double.class) {
                fileWriter.write(field.getName()+" "+field.getDouble(wrapper));
            }
            else if (type == boolean.class) {
                fileWriter.write(field.getName()+" "+field.getBoolean(wrapper));
            }
            else if (type.isArray()) {
                Class subType = type.getComponentType();
                if (subType == int.class) {
                    int[] array = (int[])field.get(wrapper);
                    for (int i=0; i<array.length; i++) {
                        fileWriter.write(field.getName()+" "+array[i]);
                    }
                }
                if (subType == long.class) {
                    long[] array = (long[])field.get(wrapper);
                    for (int i=0; i<array.length; i++) {
                        fileWriter.write(field.getName()+" "+array[i]);
                    }
                }
                else if (subType == double.class) {
                    double[] array = (double[])field.get(wrapper);
                    for (int i=0; i<array.length; i++) {
                        fileWriter.write(field.getName()+" "+array[i]);
                    }
                }
                else if (subType == boolean.class) {
                    boolean[] array = (boolean[])field.get(wrapper);
                    for (int i=0; i<array.length; i++) {
                        fileWriter.write(field.getName()+" "+array[i]);
                    }
                }
                else {
                    writeUnknownType(fileWriter, field);
                }
            }
            else {
                writeUnknownType(fileWriter, field);
            }
        }
        catch (IllegalAccessException e) {
            throw new RuntimeException("Illegal access exception trying to get "+field.getName()+" from "+wrapper, e);
        }
        fileWriter.write("\n");
    }
    
    protected void writeUnknownType(FileWriter fileWriter, Field field) {
        if (wrapper instanceof ObjectParamWrapper) {
            ((ObjectParamWrapper)wrapper).writeField(fileWriter, field);
        }
        else {
            throw new RuntimeException("don't know how to parse field "+field.getName());
        }
    }
    
    private static final long serialVersionUID = 1L;
    protected ParamBase wrapper;
    protected String fileName;
}
