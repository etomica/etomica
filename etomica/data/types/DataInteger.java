package etomica.data.types;

import etomica.data.Data;
import etomica.data.DataInfo;
import etomica.data.DataInfoFactory;
import etomica.data.types.DataArithmetic.DataInfoArithmetic;
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
public class DataInteger implements Data, java.io.Serializable {

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
     * The encapsulated <tt>int</tt> value.
     */
    public int x;
    
    public static class DataInfoInteger extends DataInfo implements DataInfoArithmetic {
        public DataInfoInteger(String label, Dimension dimension) {
            super(label, dimension);
        }
        
        public DataInfoFactory getFactory() {
            return new DataInfoIntegerFactory(this);
        }
    }
    
    public static class DataInfoIntegerFactory extends DataInfoFactory {
        protected DataInfoIntegerFactory(DataInfoInteger template) {
            super(template);
        }
        
        public DataInfo makeDataInfo() {
            return new DataInfoInteger(label, dimension);
        }
    }
}
