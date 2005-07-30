package etomica.data.types;

import etomica.Data;
import etomica.DataInfo;
import etomica.data.DataFactory;
import etomica.units.Dimension;


/**
 * Data object encapsulating a single mutable value of type integer. <br>
 * 
 * @author Andrew Schultz and David Kofke
 *  
 */

/*
 * History
 * Created on Jun 15, 2005 by kofke
 */
public class DataInteger extends Data {

    /**
     * Constructs a new instance with the given DataInfo.
     */
    public DataInteger(String label, Dimension dimension) {
        super(new DataInfo(label, dimension, getFactory()));
    }

    /**
     * Copy constructor.
     */
    public DataInteger(DataInteger data) {
        super(data);
        x = data.x;
    }
    
    /**
     * Returns a copy of this instance.  Returned object has its own instances of
     * all fields, set equal to the values of this instance's fields.
     */
    public Data makeCopy() {
        return new DataInteger(this);
    }

    /**
     * Sets the encapsulated integer to the given value.
     */
    public void E(Data y) {
        x = ((DataInteger)y).x;
    }

    /**
     * Sets the encapsulated integer to the given value.
     */
    public void E(int y) {
        x = y;
    }
    
    /**
     * Returns a string formed from the dataInfo label and the encapsulated integer.
     */
    public String toString() {
        return dataInfo.getLabel() + " " + Integer.toString(x);
    }
    
    /**
     * Returns the singleton instance of a factory that makes this data type.
     */
    public static DataFactory getFactory() {
        return FACTORY;
    }
    
    public int x;
    
    private static final Factory FACTORY = new Factory();

    private static class Factory implements DataFactory {
        
        public Data makeData(String label, Dimension dimension) {
            return new DataInteger(label, dimension);
        }
        
        public Class getDataClass() {
            return DataInteger.class;
        }
        
    }

}
