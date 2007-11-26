package etomica.graphics;

import java.util.Vector;

import javax.swing.event.TableModelEvent;
import javax.swing.event.TableModelListener;
import javax.swing.table.AbstractTableModel;

public class DeviceTableModelGeneric extends AbstractTableModel {

    private String[] columnNames;
    private Class[] columnClasses;
    private Vector data;
    private int numColumns = 0;
    private int numRows = 0;

    public DeviceTableModelGeneric(Object[][] initData, String[] colNames) {
        columnNames = colNames;

        if(initData != null) {
        	if(colNames.length != initData[0].length) {
        		System.out.println("WARNING : GENERATE AN EXCEPTION");
        	}
            numRows = initData.length;
            numColumns = initData[0].length;
        }
        else {
        	numRows = 0;
        	numColumns = colNames.length;
        }

        columnClasses = new Class[numColumns];
        for(int i = 0; i < numColumns; i++) {
        	columnClasses[i] = Number.class;
        }

        data = new Vector();
        if(initData != null) {
            for(int row = 0; row < numRows; row++) {
        	    data.add(new TableDataType(row, numColumns, initData[row]));
            }
        }
        addEmptyRow();

        addTableModelListener(new TableDataListener());
    }
    

    /**
     * Get the value in the given table cell.
     * @param row row where cell resides
     * @param col column where cell resides
     * @return the value of the given table cell wrapped as an Object
     */
    public Object getValueAt(int row, int col) {
        if (row == TableModelEvent.HEADER_ROW)
            return columnNames[col];
        else if (col == TableModelEvent.ALL_COLUMNS)
        	return null;

    	Object value = "";

    	if(data.size() >= row) {
    	    value = ((TableDataType)(data.get(row))).value[col];
    	}
    	return value;
    }

    /**
     * Set the value in the given table cell.
     * @param the value to set in the table cell wrapped as an Object
     * @param row row where cell to set resides
     * @param col column where cell to set resides
     */
    public void setValueAt(Object value, int row, int col) {
    	((TableDataType)data.get(row)).value[col] = value;
    	fireTableCellUpdated(row, col);
    }

    /**
     * Get the number of rows in the table (including blank rows).
     */
    public int getRowCount() {
    	return numRows;
    }

    /**
     * Get the number of columns in the table
     */
    public int getColumnCount() {
    	return numColumns;
    }

    /**
     * Delete an entire row from the table
     * @param row row to delete from the table.
     */
    public void deleteRow(int row) {
        data.remove(row);
    	numRows--;
        fireTableRowsDeleted(row, row);
    }

    public boolean isCellEditable(int row, int column) {
    	return true;
    }

    public String getColumnName(int column) {
    	return columnNames[column];
    }

    public void setColumnNames(String[] newColumnNames) {
        columnNames = newColumnNames;
    }

    public Class getColumnClass(int column) {
	    return columnClasses[column];
    }

    /**
     * Adds an empty row to the end of the table.
     */
    private void addEmptyRow() {
    	numRows++;
    	data.add(new TableDataType(numRows, numColumns, createBlankObject(numColumns)));
    	fireTableRowsInserted(getRowCount(), getRowCount());
    }

    private Object[] createBlankObject(int size) {
    	Object[] obj = new Object[size];
    	for(int i = 0; i < size; i++) {
    		obj[i] = "";
    	}
    	return obj;
    }

    private class TableDataType {
    	public int row;
    	public int numCols;
    	public Object[] value;
    	
    	public TableDataType(int r, int c, Object[] val) {
    		row = r;
    		numCols = c;
    		value = new Object[numCols];
    		for(int col = 0; col < numCols; col++) {
    			value[col] = val[col];
    		}
    	}
    } // end class TableDataType

    // Hmmm, no order guarenteed between this listener and other listeners...
    private class TableDataListener implements TableModelListener {

    	public TableDataListener() {
    	}

    	public void tableChanged(TableModelEvent e) {

    		Object blank = "";

            if(e.getType() == TableModelEvent.UPDATE) {
            	if(e.getColumn() == TableModelEvent.ALL_COLUMNS) {
            		;
            	}
            	else if(e.getFirstRow() == TableModelEvent.HEADER_ROW) {
            		;
            	}
            	else {
		    		if(e.getFirstRow() == numRows-1) {
		    		    if(getValueAt(e.getFirstRow(), e.getColumn()).equals(blank) == false) {
		    		    	addEmptyRow();	
		    			}
		    		}
            	}
	    	}
            else if(e.getType() == TableModelEvent.DELETE) {
            }
    	}
    } // end class TableDataListener

}
