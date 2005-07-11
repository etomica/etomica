package etomica;

import etomica.data.DataFactory;
import etomica.units.Dimension;

/**
 * Object held by a Data instance and which provides descriptive information
 * about the data encapsulated in it. Information is typically used when
 * displaying or writing the data to file.
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
     *            descriptive label for the data; may be changed after
     *            construction
     * @param dimension
     *            physical dimensions (e.g., length, force) of the data; cannot
     *            be changed after construction
     */
    public DataInfo(String label, Dimension dimension, DataFactory factory) {
        this.label = label;
        this.dimension = dimension;
        this.dataFactory = factory;
    }
    
    //remove this constructor
    public DataInfo(String label, Dimension dimension) {
        this.label = label;
        this.dimension = dimension;
        dataFactory = null;
    }

    /**
     * Copy constructor. Makes new instance with fields equal to those of the
     * given instance.
     */
    public DataInfo(DataInfo dataInfo) {
        this.label = new String(dataInfo.label);
        this.dimension = dataInfo.dimension;
        this.dataFactory = dataInfo.dataFactory.copy();
    }
    
    public DataFactory getDataFactory() {
        return dataFactory;
    }
    
    public Class getDataClass() {
        return dataFactory.getDataClass();
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
     * Sets the descriptive label for the data.
     * 
     * @param label
     *            new label
     */
    public void setLabel(String label) {
        this.label = label;
    }
    
    public String toString() {
        return label + " (" + dimension.toString() + ")";
    }

    private String label;
    private final Dimension dimension;
    private final DataFactory dataFactory;
}
