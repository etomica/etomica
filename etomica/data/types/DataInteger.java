package etomica.data.types;

import etomica.data.Data;
import etomica.data.DataInfo;
import etomica.data.DataInfoFactory;
import etomica.data.DataTag;
import etomica.data.types.DataArithmetic.DataInfoArithmetic;
import etomica.units.Dimension;


/**
 * Data object encapsulating a single mutable value of type <tt>int</tt>.
 * 
 * @author Andrew Schultz and David Kofke
 *  
 */
public class DataInteger implements Data, java.io.Serializable {

    /**
     * Constructs a new instance with the given DataInfo.
     */
    public DataInteger() {
        super();
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
    
    private static final long serialVersionUID = 1L;
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

        public Data makeData() {
            return new DataInteger();
        }

        private static final long serialVersionUID = 1L;
    }
    
    public static class DataInfoIntegerFactory extends DataInfoFactory {
        protected DataInfoIntegerFactory(DataInfoInteger template) {
            super(template);
        }
        
        public DataInfo makeDataInfo() {
            DataInfoInteger dataInfo = new DataInfoInteger(label, dimension);
            DataTag[] tagArray = new DataTag[tags.size()];
            dataInfo.addTags((DataTag[])tags.toArray(tagArray));
            return dataInfo;
        }

        private static final long serialVersionUID = 1L;
    }
}
