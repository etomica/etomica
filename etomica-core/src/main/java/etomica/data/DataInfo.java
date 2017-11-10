/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.data;

import java.util.ArrayList;

import etomica.units.dimensions.Dimension;
import etomica.util.Debug;

/**
 * Object associated with a Data instance and which provides descriptive 
 * information about the data encapsulated in it. Information is typically used
 * when displaying or writing the data to file, or making new Data instances 
 * that can work with the one holding this DataInfo.
 * <p>
 * A DataInfo instance is completely immutable.
 * 
 * @author Andrew Schultz and David Kofke
 * 
 * @see IData
 */
public abstract class DataInfo implements java.io.Serializable, IDataInfo {

    /**
     * Constructs new instance with descriptive label and dimension.
     * 
     * @param label
     *            descriptive label for the data
     * @param dimension
     *            physical dimensions (e.g., length, force) of the data
     */
    protected DataInfo(String label, Dimension dimension) {
        this.label = label;
        this.dimension = dimension;
        tags = new ArrayList();
        tagArray = new DataTag[0];
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

    /**
     * Adds the tag
     */
    public void addTag(DataTag newTag) {
        if (Debug.ON && newTag.getClass().isArray()) {
            System.err.println("You probably wanted addTags");
        }
        tags.add(newTag);
    }
    
    /**
     * Convenience method to add all of the given tags.
     */
    public void addTags(DataTag[] newTags) {
        for (int i=0; i<newTags.length; i++) {
            tags.add(newTags[i]);
        }
    }
    
    /**
     * Removes all tags.
     */
    public void clearTags() {
        tags.clear();
    }
    
    /**
     * Returns the array of tags held by this DataInfo instance.  The same
     * instance of the array can be returned each time the method is called,
     * so the array elements should not be modified.
     */
    public DataTag[] getTags() {
        if (tagArray.length != tags.size()) {
            tagArray = new DataTag[tags.size()];
        }
        tagArray = (DataTag[])tags.toArray(tagArray);
        return tagArray;
    }
    
    private final String label;
    private final Dimension dimension;
    private final ArrayList tags;
    private DataTag[] tagArray;
}
