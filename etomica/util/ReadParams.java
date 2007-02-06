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
    
    public ReadParams(String inputFileName, ParamBase parameterWrapper) {
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
        fields = wrapper.getClass().getFields();
    }

    /**
     * Reads each line from inputFileName and attempts to match it with a field
     * from parameterWrapper and sets the field to the value from the file.
     * This routine handles boolean, int, long, double, String and arrays of 
     * boolean, int long and double.
     */
    public boolean readParameters() {
        if (fileName == null) {
            throw new IllegalStateException("you must set the file name to read");
        }
        boolean success = true;
        FileReader fileReader;
        try {
            fileReader = new FileReader(fileName);
        }catch(IOException e) {
            if (firstException == null) {
                firstException = e;
            }
            throw new RuntimeException("Cannot open "+fileName, e);
        }
        BufferedReader bufReader = new BufferedReader(fileReader);
        boolean classRead = false;
        
        while (true) {
            String line;
            try {
                line = bufReader.readLine();
            }catch(IOException e) {
                if (firstException == null) {
                    firstException = e;
                }
                try {
                    bufReader.close();
                }
                catch (IOException e2) {
                    // hmm.  couldn't close either.  something's b0rked
                }
                throw new RuntimeException("Cannot read "+fileName, e);
            }
            // bail on EOF
            if (line == null) break;
            if (!classRead && line.matches("^# Class [^ ].*$")) {
                // check to see if the comment line indicates the wrapper class
                classRead = true;
                int i = new String("# Class ").length();
                String className = line.substring(i).trim();
                Class wrapperClass = getWrapperClass(className);
                if (wrapperClass != null) {
                    if (wrapper == null) {
                        wrapper = makeWrapperObject(wrapperClass);
                        if (wrapper == null) {
                            System.err.println("Could not create an instance of class "+wrapperClass.getName());
                            return false;
                        }
                        setParameterWrapper(wrapper);
                    }

                    if (!wrapperClass.isAssignableFrom(wrapper.getClass())) {
                        System.err.println("Input file requires an object of class "+wrapperClass.getName()+" but you gave me "+wrapper.getClass());
                        return false;
                    }
                    continue;
                }
                // first line apparently wasn't the class name
                if (wrapper == null) {
                    throw new IllegalStateException("You must either provide the wrapper Object or have the wrapper Class name at the top of the input file");
                }
                // forget the exception we got because we couldn't find the wrapper class
                // it's optional, so long as a wrapper object is provided
                firstException = null;
            }
            else if (line.matches("^ *#.*")) {
                // skip comments
                continue;
            }
            int i = line.indexOf(' ');
            String token = line.substring(0,i).trim();
            String value = line.substring(i).trim();
            if (value.length() == 1) {
                if (firstException == null) {
                    firstException = new RuntimeException("bogus line encountered in "+fileName);
                }
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
                if (firstException == null) {
                    firstException = new RuntimeException("don't know what to do with line: "+line);
                }
                System.err.println("don't know what to do with line:");
                System.err.println(line);
                success = false;
            }
        }
        try {
            bufReader.close();
        }
        catch (IOException e) {
        }
        return success;
    }
    
    protected Class getWrapperClass(String line) {
        Class wrapperClass = null;
        try {
            wrapperClass = Class.forName(line);
        }
        catch (ClassNotFoundException e) {
            if (firstException == null) {
                firstException = e;
            }
            // if the declared class doesn't exist, we can't parse
            // the file
            return null;
        }
        return wrapperClass;
    }
    
    public ParamBase makeWrapperObject(Class wrapperClass) {
        try {
            // make sure we can instantiate this thing
            return (ParamBase)wrapperClass.newInstance();
        }
        catch (InstantiationException e) {
            if (firstException == null) {
                firstException = e;
            }
            return null;
        }
        catch (IllegalAccessException e) {
            if (firstException == null) {
                firstException = e;
            }
            return null;
        }
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
            else if (type == float.class) {
                field.setFloat(wrapper,Float.parseFloat(value));
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
                else if (subType == float.class) {
                    float[] array = new float[strings.length-1];
                    for (int i=0; i<array.length; i++) {
                        array[i] = Float.parseFloat(strings[i+1]);
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
                    return parseUnknownType(field, value);
                }
            }
            else {
                return parseUnknownType(field, value);
            }
        }
        catch (NumberFormatException ex) {
            // catches exceptions thrown by Type.parseType for
            // int, long, float and double
            if (firstException == null) {
                firstException = new RuntimeException("error parsing "+value+" for "+field.getName(), ex);
            }
            System.out.println("error parsing "+value+" for "+field.getName());
            return false;
        }
        catch (IllegalAccessException e) {
            if (firstException == null) {
                firstException = new RuntimeException("Illegal access exception trying to set "+field.getName()+" to "+value, e);
            }
            System.err.println("Illegal access exception trying to set "+field.getName()+" to "+value);
            return false;
        }
        return true;
    }
    
    protected boolean parseUnknownType(Field field, String value) {
        if (wrapper instanceof ObjectParamWrapper) {
            try {
                ((ObjectParamWrapper)wrapper).parseField(field, value);
            }
            catch (RuntimeException ex) {
                if (firstException == null) {
                    firstException = ex;
                }
                System.err.println("couldn't parse field "+field.getName()+" with value "+value);
                return false;
            }
        }
        else {
            if (firstException == null) {
                firstException = new RuntimeException("don't know how to parse field "+field.getName());
            }
            System.err.println("don't know how to parse field "+field.getName());
            return false;
        }
        return true;
    }
    
    /**
     * Returns the first exception thrown when attempting to read the file
     */
    public Exception getFirstException() {
        return firstException;
    }

    private static final long serialVersionUID = 2L;
    protected Exception firstException;
    protected Field[] fields;
    protected ParamBase wrapper;
    protected String fileName;
}
