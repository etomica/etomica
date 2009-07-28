package etomica.data;

public interface IDataInfo {

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
     * Returns the number of numerical values held by the IData object.
     */
    public int getLength();

    /**
     * Returns a Data object appropriate for this DataInfo instance.
     */
    public IData makeData();

}