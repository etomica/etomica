package etomica.data;



/**
 * Interface for a class that performs some action in response to
 * a change in a DataSet instance.  Classes implementing this interface can
 * register with a DataSet via the addDataListener method.  Typically
 * this interface is used to cause a display element to update when the
 * data is updated.
 * 
 * @author Andrew Schultz
 */
public interface DataSetListener {

    /**
     * Method called when one or more pieces of data have changed.
     */
    public void dataChanged(DataSet dataSet);
    
    /**
     * Method called when a Data object is added.  Firing is
     * performed after the object is added to the DataSet.
     * @param newData the Data object that has been added
     */
    public void dataCountChanged(DataSet dataSet);
}
