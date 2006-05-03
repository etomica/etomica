package etomica.data;

import etomica.units.Dimension;

/**
 * Object held by a Data instance and which provides descriptive information
 * about the data encapsulated in it. Information is typically used when
 * displaying or writing the data to file, or making new Data instances that can
 * work with the one holding this DataInfo.
 * <p>
 * A DataInfo instance is completely immutable, and it is declared final in the
 * Data instance holding it.
 * 
 * @author Andrew Schultz and David Kofke
 * 
 * @see Data
 *  
 */

/*
 * History Created on Jun 15, 2005 by kofke
 */
public class DataInfo implements java.io.Serializable {

    /**
     * Constructs new instance with descriptive label and dimension.
     * 
     * @param label
     *            descriptive label for the data
     * @param dimension
     *            physical dimensions (e.g., length, force) of the data
     * @param factory
     *            a DataFactory that makes Data instances of the type holding
     *            this DataInfo. New Data instances will be independent of the
     *            one holding this, but will be structured the same way
     */
    public DataInfo(String label, Dimension dimension) {
        this.label = label;
        this.dimension = dimension;
    }

    /**
     * @return Returns the dimension given at construction.
     */
    public Dimension getDimension() {
        return dimension;
    }

    /**
     * @return Returns the descriptive label of the data, as given in
     *         constructor or at last call to setLabel.
     */
    public String getLabel() {
        return label;
    }

    /**
     * Returns "label (dimension)", where label and dimension are the values
     * held by this instance.
     */
    public String toString() {
        return label + " (" + dimension.toString() + ")";
    }

    private final String label;
    private final Dimension dimension;
    
    /**
     * Marker interface for DataInfo classes that correspond to DataArithmetic classes
     * @author Andrew Schultz
     */
    public interface DataInfoArithmetic {}
}
