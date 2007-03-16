package etomica.data;

import etomica.units.Dimension;

public interface IDataInfo {

    /**
     * @return Returns the dimension given at construction.
     */
    public Dimension getDimension();

    /**
     * @return Returns the descriptive label of the data, as given in
     *         constructor or at last call to setLabel.
     */
    public String getLabel();

    /**
     * Returns "label (dimension)", where label and dimension are the values
     * held by this instance.
     */
    public String toString();

    /**
     * Adds the tag to this IDataInfo object.  You should only call this method
     * if you created this IDataInfo (You should not call this for IDataInfo
     * you receive.  Make a new instance with the DataInfoFactory and then add
     * your tag to that).
     */
    public void addTag(DataTag newTag);

    /**
     * Convenience method to add all of the given tags.  You should only call
     * this method if you created this IDataInfo (You should not call this for
     * IDataInfo you receive.  Make a new instance with the DataInfoFactory and
     * then add your tag to that).
     */
    public void addTags(DataTag[] newTags);

    /**
     * Removes all tags.  You should only call this method if you created this
     * IDataInfo.
     */
    public void clearTags();

    /**
     * Returns the array of tags held by this DataInfo instance.  The same
     * instance of the array can be returned each time the method is called,
     * so the array elements should not be modified.
     */
    public DataTag[] getTags();

    /**
     * Returns a mutable factory that can make copies of this instance of
     * DataInfo.
     */
    public DataInfoFactory getFactory();

    /**
     * Returns a Data object appropriate for this DataInfo instance.
     */
    public Data makeData();

}