package etomica.data;

import java.util.ArrayList;

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
public abstract class DataInfo implements java.io.Serializable {

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
    protected DataInfo(String label, Dimension dimension) {
        this.label = label;
        this.dimension = dimension;
        tags = new ArrayList();
        tagArray = new Object[0];
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
    
    public void addTag(Object newTag) {
        tags.add(newTag);
    }
    
    public void addTags(Object[] newTags) {
        for (int i=0; i<newTags.length; i++) {
            tags.add(newTags[i]);
        }
    }
    
    public void clearTags() {
        tags.clear();
    }
    
    public Object[] getTags() {
        tagArray = tags.toArray(tagArray);
        return tagArray;
    }
    
    public abstract DataInfoFactory getFactory();

    private final String label;
    private final Dimension dimension;
    private ArrayList tags;
    private Object[] tagArray;
}
