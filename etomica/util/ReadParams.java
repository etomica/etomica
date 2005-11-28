package etomica.util;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.lang.reflect.Field;

/**
 * This class reads input parameters from a file and assigns the values to the
 * fields of an object.
 */
public class ReadParams implements java.io.Serializable {

    public ReadParams() {}

    public ReadParams(String inputFileName) {
        this();
        setInputFileName(inputFileName);
    }
    
    public ReadParams(String inputFileName, Object parameterWrapper) {
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
    public Object getParameterWrapper() {
        return wrapper;
    }

    /**
     * Sets the parameterWrapper
     */
    public void setParameterWrapper(Object newParameterWrapper) {
        wrapper = newParameterWrapper;
        fields = wrapper.getClass().getFields();
    }

    /**
     * Reads each line from inputFileName and attempts to match it with a field
     * from parameterWrapper and sets the field to the value from the file.
     * This routine handles boolean, int, long, double, String and arrays of 
     * boolean, int long and double.
     */
    public boolean readParameters() {
        boolean success = true;
        FileReader fileReader;
        try {
            fileReader = new FileReader(fileName);
        }catch(IOException e) {
            throw new RuntimeException("Cannot open "+fileName+", caught IOException: " + e.getMessage());
        }
        BufferedReader bufReader = new BufferedReader(fileReader);
        
        while (true) {
            String line;
            try {
                line = bufReader.readLine();
            }catch(IOException e) {
                throw new RuntimeException("Cannot read "+fileName+", caught IOException: " + e.getMessage());
            }
            // bail on EOF
            if (line == null) break;
            // skip comments
            if (line.matches("^ *#.*$")) continue;
            int i = line.indexOf(' ');
            String token = line.substring(0,i).trim();
            String value = line.substring(i).trim();
            if (value.length() == 1) {
                System.err.println("bogus line encountered in "+fileName);
                System.err.println(line);
                success = false;
                continue;
            }
            boolean foundField = false;
            for (int j=0; j<fields.length; j++) {
                if (token.equals(fields[j].getName())) {
                    success = setValue(fields[j],value) && success;
                    foundField = true;
                    break;
                }
            }
            if (!foundField) {
                System.err.println("don't know what to do with line:");
                System.err.println(line);
                success = false;
            }
        }
        return success;
    }
    
    protected boolean setValue(Field field, String value) {
        Class type = field.getType();
        try {
            if (type == String.class) {
                field.set(wrapper,value);
            }
            else if (type == int.class) {
                field.setInt(wrapper,Integer.parseInt(value));
            }
            else if (type == long.class) {
                field.setLong(wrapper,Long.parseLong(value));
            }
            else if (type == double.class) {
                field.setDouble(wrapper,Double.parseDouble(value));
            }
            else if (type == boolean.class) {
                field.setBoolean(wrapper,Boolean.valueOf(value).booleanValue());
            }
            else if (type.isArray()) {
                Class subType = type.getComponentType();
                String[] strings = value.split(" +");
                if (subType == int.class) {
                    int[] array = new int[strings.length-1];
                    for (int i=0; i<array.length; i++) {
                        array[i] = Integer.parseInt(strings[i+1]);
                    }
                    field.set(wrapper,array);
                }
                if (subType == long.class) {
                    long[] array = new long[strings.length-1];
                    for (int i=0; i<array.length; i++) {
                        array[i] = Long.parseLong(strings[i+1]);
                    }
                    field.set(wrapper,array);
                }
                else if (subType == double.class) {
                    double[] array = new double[strings.length-1];
                    for (int i=0; i<array.length; i++) {
                        array[i] = Double.parseDouble(strings[i+1]);
                    }
                    field.set(wrapper,array);
                }
                else if (subType == boolean.class) {
                    boolean[] array = new boolean[strings.length-1];
                    for (int i=0; i<array.length; i++) {
                        array[i] = Boolean.valueOf(strings[i+1]).booleanValue();
                    }
                    field.set(wrapper,array);
                }
                else {
                    System.err.println("unknown array type "+subType+" for "+field.getName());
                    return false;
                }
            }
            else {
                System.err.println("unknown type "+type+" for "+field.getName());
                return false;
            }
        }
        catch (IllegalAccessException e) {
            System.err.println("Illegal access exception trying to set "+field.getName()+" to "+value);
            return false;
        }
        return true;
    }

    protected Field[] fields;
    protected Object wrapper;
    protected String fileName;
}
