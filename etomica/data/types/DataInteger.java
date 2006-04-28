package etomica.data.types;

import java.io.Serializable;

import etomica.data.Data;
import etomica.data.DataFactory;
import etomica.data.DataInfo;
import etomica.units.Dimension;


/**
 * Data object encapsulating a single mutable value of type <tt>int</tt>.
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
    public DataInteger() {
        super();
    }

    /**
     * Copy constructor.
     */
    public DataInteger(DataInteger data) {
        super();
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
        return Integer.toString(x);
    }
    
    /**
     * Returns a (singleton) DataFactory that makes DataInteger instances.
     */
    public static DataFactory getFactory() {
        return FACTORY;
    }
    
    /**
     * The encapsulated <tt>int</tt> value.
     */
    public int x;
    
    private static final Factory FACTORY = new Factory();

    private static class Factory implements DataFactory, Serializable {
        
        public Data makeData() {
            return new DataInteger();
        }
        
        public Class getDataClass() {
            return DataInteger.class;
        }
        
    }

}
