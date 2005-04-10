package etomica.data;


/**
 * Interface for a class that performs some action in response to
 * a change in a DataTable instance.  Classes implementing this interface can
 * register with a DataTable via the addTableListener method.  Typically
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
public interface DataTableListener {

    /**
     * Method called when one or more entries in the table have changed.
     */
    public void tableDataChanged();
    
    /**
     * Method called when a table column is added.  Firing is
     * performed after the column is added to the table.
     * @param newColumn the column that has been added
     */
    public void tableColumnAdded(DataBin newColumn);
    
    /**
     * Method called when a table column is removed.  Firing is
     * performed after column is removed from table.
     * @param index index of the column before it was removed
     * @param removedColumn handle to the removed column
     */
    public void tableColumnRemoved(int index, DataBin removedColumn);

    /**
     * Method called to indicate that the number of rows in the
     * table has changed.  Method is called after the change takes place.
     * @param oldCount the old number of rows
     * @param newCount the new number of rows
     */
    public void tableRowCountChanged(int oldCount, int newCount);
}
