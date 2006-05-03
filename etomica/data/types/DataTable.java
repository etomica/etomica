package etomica.data.types;

import java.io.Serializable;

import etomica.data.Data;
import etomica.data.DataInfo;
import etomica.units.Null;
import etomica.util.Function;

/**
 * Data object that holds <tt>double[]</tt> arrays as if they are columns in a
 * table. The inner class <tt>Column</tt> wraps the array along with a
 * descriptive String heading and a Dimension instance that indicates the
 * physical dimensions of the data in the column.  Descriptive immutable String 
 * labels can also be associated with the rows (row headers).
 * be 
 * <p>
 * The number of columns is set at construction and cannot be changed.
 * <p>
 * All columns are of equal length (i.e., the length of the <tt>double[]</tt>
 * arrays is the same for all columns). The common length of the columns is
 * set at construction and cannot be changed.
 * <p>
 * The label in the DataInfo for this object should be descriptive of the table
 * as a whole; the dimension in the DataInfo will be Dimension.MIXED if the dimensions
 * of the data in the columns are not all the same, otherwise it will be the common
 * dimension of the data in the columns.
 * 
 * @author David Kofke and Andrew Schultz
 *  
 */

/*
 * History Created on Jul 28, 2005 by kofke
 */
public class DataTable extends DataGroup implements DataArithmetic, Serializable {

    /**
     * Creates a new table with a specified number of columns all of a given
     * length. The assigned heading for each column is its ordinal (starting
     * from 0), e.g., the heading of the first column is '0'.  All columns are
     * assigned the given Dimension, which is also the Dimension of the DataTable
     * (as kept in DataInfo).  No row headers are defined.
     * 
     * @param label
     *            a descriptive label for the table
     * @param dimension
     *            indicates the common Dimension for all columns
     * @param nColumns
     *            the number of columns to make
     * @param nRows
     *            the length of each and every column
     */
    public DataTable(int nColumns, int nRows) {
        super(makeData(nColumns, nRows));
    }

    protected static DataDoubleArray[] makeData(int nColumns, int nRows) {
        DataDoubleArray[] data = new DataDoubleArray[nColumns];
        for (int i=0; i<nColumns; i++) {
            data[i] = new DataDoubleArray(nRows);
        }
        return data;
    }
    
    public DataTable(DataDoubleArray[] columns) {
        super(columns);
        for (int i=1; i<columns.length; i++) {
            if (columns[i].getLength() != columns[0].getLength()) {
                throw new IllegalArgumentException("All columns must have the same length");
            }
        }
    }
    
    /**
     * Copy constructor. Makes a new DataTable with new Column instances having
     * the same data values and row headings as in the given DataTable.
     */
    public DataTable(DataTable dataTable) {
        this(dataTable.getNData(), dataTable.getNRows());
        E(dataTable);
    }

    /**
     * Makes a copy of this DataTable using the copy constructor.
     */
    public Data makeCopy() {
        return new DataTable(this);
    }
    
    /**
     * Copies all data values from the given table to this one.
     * 
     * @throws ClassCastException
     *             if the argument is not an instance of DataTable
     */
    public void E(Data otherData) {
        DataTable table = (DataTable) otherData;
        for (int i = 0; i < data.length; i++) {
            ((DataDoubleArray)data[i]).E((DataDoubleArray)table.getData(i));
        }
    }

    /**
     * Returns the number of rows in each and every column of the table.
     * 
     * @throws ArrayIndexOutOfBoundsException
     *             if the DataTable has zero columns
     */
    public int getNRows() {
        return ((DataDoubleArray)data[0]).getLength();
    }
    
    /**
     * Plus-equals (+=) operation.
     */
    public void PE(DataArithmetic y) {
        DataTable table = (DataTable) y;
        for (int i = 0; i < data.length; i++) {
            ((DataDoubleArray)data[i]).PE((DataDoubleArray)table.getData(i));
        }
    }

    /**
     * Minus-equals (-=) operation.
     */
    public void ME(DataArithmetic y) {
        DataTable table = (DataTable) y;
        for (int i = 0; i < data.length; i++) {
            ((DataDoubleArray)data[i]).ME((DataDoubleArray)table.getData(i));
        }
    }

    /**
     * Times-equals (*=) operation.
     */
    public void TE(DataArithmetic y) {
        DataTable table = (DataTable) y;
        for (int i = 0; i < data.length; i++) {
            ((DataDoubleArray)data[i]).TE((DataDoubleArray)table.getData(i));
        }
    }

    /**
     * Divide-equals (/=) operation.
     */
    public void DE(DataArithmetic y) {
        DataTable table = (DataTable) y;
        for (int i = 0; i < data.length; i++) {
            ((DataDoubleArray)data[i]).DE((DataDoubleArray)table.getData(i));
        }
    }

    /**
     * Equals (=) operation, sets all values in data equal to the given value.
     */
    public void E(double y) {
        for (int i = 0; i < data.length; i++) {
            ((DataDoubleArray)data[i]).E(y);
        }
    }

    /**
     * Plus-equals (+=) operation, adding given value to all values in data.
     */
    public void PE(double y) {
        for (int i = 0; i < data.length; i++) {
            ((DataDoubleArray)data[i]).PE(y);
        }
    }

    /**
     * Times-equals (*=) operation, multiplying all values in data by given
     * value.
     */
    public void TE(double y) {
        for (int i = 0; i < data.length; i++) {
            ((DataDoubleArray)data[i]).TE(y);
        }
    }

    /**
     * Maps the function on all data values, replace each with the value given
     * by the function applied to it.
     */
    public void map(Function function) {
        for (int i = 0; i < data.length; i++) {
            ((DataDoubleArray)data[i]).map(function);
        }
    }

    /**
     * Returns the number of values held by the data instance.
     */
    public int getLength() {
        if(data.length == 0) {
            return 0;
        }
        return getNRows() * data.length;
    }

    /**
     * Returns the i-th data value, where the values are numbered by counting down the first column,
     * then the next, and so on. 
     * <p>
     * Specifically, returns column[i/nRows].data[i % nRows].
     * 
     * @param i
     *            index of the desired data value
     * @return the data value
     * @throws IllegalArgumentException
     *             if i >= getLength()
     */
    public double getValue(int i) {
        if(i < 0 || i >= getLength()) {
            throw new IllegalArgumentException("Illegal index for get value. Maximum is "+getLength());
        }
        int nRows = getNRows();
        return ((DataDoubleArray)data[i/nRows]).getValue(i % nRows);
    }

    /**
     * Returns a new array formed from the values held by the data instance.
     * 
     * @throws IllegalArgumentException
     *             if the length of the array is not equal to getLength
     */
    public void assignTo(double[] array) {
        if(array.length != getLength()) {
            throw new IllegalArgumentException("Illegal array size; length must be: "+getLength());
        }
        if(data.length == 0) {
            return;
        }
        int nRows = getNRows();
        for(int i=0; i<data.length; i++) {
            System.arraycopy(((DataDoubleArray)data[i]).getData(), 0, array, i*nRows, nRows);
        }
    }

    /**
     * Returns true if any data value is true for Double.isNaN
     */
    public boolean isNaN() {
        for (int i = 0; i < data.length; i++) {
            if (((DataDoubleArray)data[i]).isNaN()) {
                return true;
            }
        }
        return false;
    }

    public static class DataInfoTable extends DataInfoGroup {
        public DataInfoTable(String label, DataInfo[] columnInfo, int nRows, String[] rowHeaders) {
            super(label, Null.DIMENSION, columnInfo);
            this.nRows = nRows;
            if (rowHeaders != null) {
                this.rowHeaders = (String[])rowHeaders.clone();
            }
            else {
                this.rowHeaders = null;
            }
        }
        
        public int getNRows() {
            return nRows;
        }
        
        public String getRowHeader(int i) {
            return rowHeaders == null ? "" : rowHeaders[i];
        }
        
        private final String[] rowHeaders;
        private final int nRows;
    }
    
}

