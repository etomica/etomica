package etomica.util;

import java.io.FileWriter;
import java.lang.reflect.Field;

/**
 * Optional interface for parameter wrappers that allows more flexibility
 * in what can be read / written to the file.  {Read,Write}Params classes
 * can handle booleans, ints, longs, floats, doubles, arrays of those types 
 * and Strings.  Other types are passed off to the wrapper object itself if
 * it is of type ParamWrapper.
 */
public interface ObjectParamWrapper {

    /**
     * Writes the given field using the given fileWriter.  The field is written
     * on a single line in the format "fieldName fieldValue". "fieldValue" 
     * corresponds to the object's current value for the given field, may 
     * contain spaces and must be parseable by the parseField method. This 
     * method should throw an exception if it is unable to write the field.
     */
    public void writeField(FileWriter fileWriter, Field field);
    
    /**
     * Parses the given field using the given fileWriter and assigns the 
     * resulting object to the field variable.  The value String should be a 
     * single line in the format written by the writeField method. This method 
     * should throw an exception if it is unable to parse the field value and
     * assign it to the field variable.
     */
    public void parseField(Field field, String value); 
}
