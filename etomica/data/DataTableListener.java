package etomica.data;



/**
 * Interface for a class that performs some action in response to
 * a change in a DataSinkTable instance.  Classes implementing this interface can
 * register with a DataSinkTable via the addTableListener method.  Typically
 * this interface is used to cause a display element to update when a
 * table entry changes.
 * 
 * @author David Kofke
 *
 */

/*
 * History
 * Created on Apr 9, 2005 by kofke
 */
public interface DataTableListener extends DataSetListener {

    /**
     * Method called to indicate that the number of rows in the
     * table has changed.  Method is called after the change takes place.
     * @param oldCount the old number of rows
     * @param newCount the new number of rows
     */
    public void tableRowCountChanged(DataSinkTable table);
}
