package etomica.data;

import etomica.data.types.DataArithmetic;
import etomica.utility.Arrays;

/**
 * Receives data from multiple streams and organizes it in the form of a table,
 * with each stream of data forming a column of the table. Fires event when data
 * in table is changed; can be configured to fire when any column of data
 * changes, or to perform action only when all columns are changed.
 * 
 * @author David Kofke
 *  
 */

/*
 * History Created on Apr 9, 2005 by kofke
 */
public class DataTable implements DataBinManager {

    public DataTable() {
    }
    
    public double getValue(int row, int column) {
        return ((DataArithmetic)columns[column].getData()).getValue(row);
    }

    public int getColumnCount() {
        return columns.length;
    }

    /**
     * Returns the length of the longest column.
     */
    public int getRowCount() {
        return rowCount;
    }

    public DataBin getColumn(int i) {
        return columns[i];
    }

    /**
     * Descriptive string for each column.  Returns an
     * array of length equal to the number of columns in
     * the table.  Each entry is obtained from the getLabel
     * method for the column. Never returns null.
     */
    public String[] getColumnHeadings() {
        //need to refresh with each call because columnHeadings isn't 
        //notified when setLabel of DataBin is called
        for (int i = 0; i < columns.length; i++) {
            columnHeadings[i] = columns[i].getDataInfo().getLabel();
        }
        return columnHeadings;
    }

    public DataBin makeColumn() {
        DataBin newColumn = new DataBin(this);
        columns = (DataBin[]) Arrays.addObject(columns, newColumn);
        update();
        fireColumnAddedEvent(newColumn);
        return newColumn;
    }

    /**
     * Removes column from table.  Gives a warning if specified column is not
     * part of the table.
     * @param oldColumn column to be removed
     */
    public void removeColumn(DataBin oldColumn) {
        int i = getColumnIndex(oldColumn);
        if(i == -1) return;
        columns = (DataBin[]) Arrays.removeObject(columns, oldColumn);
        update();
        fireColumnRemovedEvent(i, oldColumn);
    }

    /**
     * Removes all columns from table, one at a time, beginning with the last one.
     *
     */
    public void removeAllColumns() {
        for(int i=columns.length-1; i>=0; i--) {
            removeColumn(columns[i]);
        }
    }

    /**
     * Receives notification that data has changed for the given bin.
     * Checks whether this has caused a changed to the row count, and 
     * fires event if it has.  Fires event to indicate data change 
     * according to protocol indicated by updatingWithAnyChange flag.
     */
    public void dataChangeNotify(DataBin bin) {
        //check for change to number of rows (length of longest column)
        int binLength = ((DataArithmetic)bin.getData()).getLength();
        if(bin == longestColumn && binLength < rowCount) {
            int oldRowCount = rowCount;
            updateRowCount();
            if(oldRowCount != rowCount) {
                fireRowCountChangedEvent(oldRowCount, rowCount);
            }
        } else if(binLength > rowCount) {
            int oldRowCount = rowCount;
            rowCount = binLength;
            longestColumn = bin;
            fireRowCountChangedEvent(oldRowCount, rowCount);
        }
        //fire change event if indicated
        if (updatingOnAnyChange) {//fire event if any column changes
            bin.dataChanged = false;
            fireDataChangedEvent();
        } else if (allColumnsChanged()) {//fire event if all columns changed
            for (int i = 0; i < columns.length; i++) {
                columns[i].dataChanged = false;//reset dataChanged flags
            }
            fireDataChangedEvent();
        }
    }
    
//    public void parameterChangeNotify(DataBin bin) {
//        int i = getColumnIndex(bin);
//        if(i == -1) return;
//        if(units[i].
//    }

    //used by dataChangeNotify, returns true if all columns have dataChanged flag true
    private boolean allColumnsChanged() {
        for (int i = 0; i < columns.length; i++) {
            if (!columns[i].dataChanged)
                return false;
        }
        return true;
    }

    /**
     * @return flag defining protocol for firing table events that
     * notify of column data changes.
     */
    public boolean isUpdatingOnAnyChange() {
        return updatingOnAnyChange;
    }

    /**
     * Describes protocol for firing table events that notify of column data
     * changes. If true, event is fired when any column notifies that it has
     * changed; if false, event is fired only when all columns have
     * notified that they have changed (since last time event was fired).
     */
    public void setUpdatingOnAnyChange(boolean updatingWithAnyChange) {
        this.updatingOnAnyChange = updatingWithAnyChange;
    }

    protected void fireDataChangedEvent() {
        for (int i = 0; i < listeners.length; i++) {
            listeners[i].tableDataChanged(this);
        }
    }

    protected void fireColumnAddedEvent(DataBin newColumn) {
        for (int i = 0; i < listeners.length; i++) {
            listeners[i].tableColumnAdded(this, newColumn);
        }
    }
    protected void fireColumnRemovedEvent(int index, DataBin removedColumn) {
        for (int i = 0; i < listeners.length; i++) {
            listeners[i].tableColumnRemoved(this, index, removedColumn);
        }
    }
    protected void fireRowCountChangedEvent(int oldCount, int newCount) {
        for (int i = 0; i < listeners.length; i++) {
            listeners[i].tableRowCountChanged(this, oldCount, newCount);
        }
    }
    
    public void addTableListener(DataTableListener newListener) {
        listeners = (DataTableListener[]) Arrays.addObject(listeners,
                newListener);
    }

    public void removeTableListener(DataTableListener oldListener) {
        listeners = (DataTableListener[]) Arrays.removeObject(listeners,
                oldListener);
    }

    /**
     * Updates table to account for addition or removal of columns.
     */
    private void update() {
        columnHeadings = new String[columns.length];
        updateRowCount();
    }

    private void updateRowCount() {
        int rowCount = 0;
        for (int i = 0; i < columns.length; i++) {
            int n = ((DataArithmetic)columns[i].getData()).getLength();
            if (n > rowCount) {
                rowCount = n;
                longestColumn = columns[i];
            }
        }
    }
    
    /**
     * Returns the index of the given column, or -1 if column is not
     * in table.
     */
    private int getColumnIndex(DataBin column) {
        for(int i=0; i<columns.length; i++) {
            if(columns[i] == column) return i;
        }
        System.err.println("Warning: attempted to access a column not contained in table.");
        return -1;
    }

    private DataBin[] columns = new DataBin[0];
    private String[] columnHeadings = new String[0];
    private DataTableListener[] listeners = new DataTableListener[0];
    private boolean updatingOnAnyChange;
    private int rowCount;
    private DataBin longestColumn = null;
}