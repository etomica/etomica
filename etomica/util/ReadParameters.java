package etomica.util;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.lang.reflect.Field;

/**
 * This class reads input parameters from a file and assigns the values to the
 * fields of an object.
 */
public class ReadParameters implements java.io.Serializable {

    public ReadParameters() {}

    public ReadParameters(String inputFileName) {
        this();
        setInputFileName(inputFileName);
    }
    
    public ReadParameters(String inputFileName, ParameterBase parameterWrapper) {
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
        fields = wrapper.getClass().getFields();
    }

    /**
     * Reads each line from inputFileName and attempts to match it with a field
     * from parameterWrapper and sets the field to the value from the file.
     * This routine handles boolean, int, long, double, String and arrays of 
     * boolean, int long and double.
     */
    public void readParameters() {
        if (fileName == null) {
            throw new IllegalStateException("you must set the file name to read");
        }
        FileReader fileReader;
        try {
            fileReader = new FileReader(fileName);
        }catch(IOException e) {
            throw new RuntimeException("Cannot open "+fileName, e);
        }
        BufferedReader bufReader = new BufferedReader(fileReader);
        boolean classRead = false;
        
        while (true) {
            String line;
            try {
                line = bufReader.readLine();
            }catch(IOException e) {
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
                            throw new RuntimeException("Could not create an instance of class "+wrapperClass.getName());
                        }
                        setParameterWrapper(wrapper);
                    }

                    if (!wrapperClass.isAssignableFrom(wrapper.getClass())) {
                        throw new RuntimeException("Input file requires an object of class "+wrapperClass.getName()+" but you gave me "+wrapper.getClass());
                    }
                    continue;
                }
                // first line apparently wasn't the class name
                if (wrapper == null) {
                    throw new IllegalStateException("You must either provide the wrapper Object or have the wrapper Class name at the top of the input file");
                }
            }
            else if (line.matches("^ *#.*") || line.length() == 0) {
                // skip comments and empty lines
                continue;
            }
            int i = line.indexOf(' ');
            String token = line.substring(0,i).trim();
            String value = line.substring(i).trim();
            if (value.length() == 0) {
                throw new RuntimeException("bogus line encountered in "+fileName+" "+line);
            }
            boolean foundField = false;
            for (int j=0; j<fields.length; j++) {
                if (token.equals(fields[j].getName())) {
                    wrapper.setValue(fields[j],value);
                    foundField = true;
                    break;
                }
            }
            if (!foundField) {
                throw new RuntimeException("don't know what to do with line: "+line);
            }
        }
        try {
            bufReader.close();
        }
        catch (IOException e) {
        }
    }
    
    protected Class getWrapperClass(String line) {
        Class wrapperClass = null;
        try {
            wrapperClass = Class.forName(line);
        }
        catch (ClassNotFoundException e) {
            throw new RuntimeException(e);
        }
        return wrapperClass;
    }
    
    public ParameterBase makeWrapperObject(Class wrapperClass) {
        try {
            // make sure we can instantiate this thing
            return (ParameterBase)wrapperClass.newInstance();
        }
        catch (InstantiationException e) {
            throw new RuntimeException(e);
        }
        catch (IllegalAccessException e) {
            throw new RuntimeException(e);
        }
    }
    
    private static final long serialVersionUID = 2L;
    protected Field[] fields;
    protected ParameterBase wrapper;
    protected String fileName;
}
