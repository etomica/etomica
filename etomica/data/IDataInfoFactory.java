package etomica.data;

import etomica.units.Dimension;

public interface IDataInfoFactory {

    /**
     * Creates a new DataInfo object using the information held by this factory.
     */
    public IDataInfo makeDataInfo();

    /**
     * Sets the label
     */
    public void setLabel(String newLabel);

    /**
     * Returns the label
     */
    public String getLabel();

    /**
     * Sets the dimension
     */
    public void setDimension(Dimension newDimension);

    /**
     * Returns the dimension
     */
    public Dimension getDimension();

    /**
     * Returns the factory's tags as an array.  The array and its elements
     * should not be modified.  If modifications are needed, setTags should
     * be used.
     */
    public DataTag[] getTags();

    /**
     * Sets the factory's tags to those in the given array.
     */
    public void setTags(DataTag[] newTags);

}