/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.data.types;

import java.io.Serializable;

import etomica.data.DataTag;
import etomica.data.IData;
import etomica.data.IEtomicaDataInfo;
import etomica.data.IEtomicaDataInfoFactory;
import etomica.data.types.DataDoubleArray.DataInfoDoubleArray;
import etomica.data.types.DataDoubleArray.DataInfoDoubleArrayFactory;
import etomica.units.dimensions.Null;
import etomica.util.Arrays;

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
public class DataTable extends DataGroup implements IData, Serializable {

    /**
     * Creates a new table with a specified number of columns all of a given
     * length. The assigned heading for each column is its ordinal (starting
     * from 0), e.g., the heading of the first column is '0'.
     * 
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
     * Returns the number of rows in each and every column of the table.
     * 
     * @throws ArrayIndexOutOfBoundsException
     *             if the DataTable has zero columns
     */
    public int getNRows() {
        return data[0].getLength();
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
        return data[i/nRows].getValue(i % nRows);
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

    private static final long serialVersionUID = 1L;

    public static class DataInfoTable extends DataInfoGroup {
        public DataInfoTable(String label, IEtomicaDataInfo[] columnInfo, int nRows, String[] rowHeaders) {
            super(label, Null.DIMENSION, columnInfo);
            this.nRows = nRows;
            if (rowHeaders != null) {
                this.rowHeaders = rowHeaders.clone();
            }
            else {
                this.rowHeaders = null;
            }
        }
        
        public IEtomicaDataInfoFactory getFactory() {
            return new DataInfoTableFactory(this);
        }
        
        /**
         * Returns the number of rows in the table
         */
        public int getNRows() {
            return nRows;
        }

        /**
         * Returns the row header for row i or an empty string if there are no
         * row headers.
         */
        public String getRowHeader(int i) {
            return rowHeaders == null ? "" : rowHeaders[i];
        }

        /**
         * Returns true if row headers exist.
         */
        public boolean hasRowHeaders() {
            return rowHeaders != null;
        }
        
        public IData makeData() {
            return new DataTable(subDataInfo.length, nRows);
        }

        private static final long serialVersionUID = 1L;
        protected final String[] rowHeaders;
        protected final int nRows;
    }
    
    public static class DataInfoTableFactory extends DataInfoGroupFactory {
        protected DataInfoTableFactory(DataInfoTable template) {
            super(template);
            if (template.rowHeaders != null) {
                rowHeaders = template.rowHeaders.clone();
            }
            nRows = template.nRows;
            columnInfoFactories = new DataInfoDoubleArrayFactory[template.subDataInfo.length];
            for (int i=0; i<columnInfoFactories.length; i++) {
                columnInfoFactories[i] = (DataInfoDoubleArrayFactory)template.subDataInfo[i].getFactory();
            }
        }
        
        public IEtomicaDataInfo makeDataInfo() {
            DataInfoDoubleArray[] columnInfo = new DataInfoDoubleArray[columnInfoFactories.length];
            for (int i=0; i<columnInfo.length; i++) {
                columnInfo[i] = (DataInfoDoubleArray)columnInfoFactories[i].makeDataInfo();
            }
            DataInfoTable dataInfo = new DataInfoTable(label, columnInfo, nRows, rowHeaders);
            DataTag[] tagArray = new DataTag[tags.size()];
            dataInfo.addTags((DataTag[])tags.toArray(tagArray));
            return dataInfo;
        }
        
        /**
         * Sets the number of rows.  If row headers exist, the array of row 
         * will be resized to contain the given number of rows.
         */
        public void setNRows(int newNRows) {
            if (nRows == newNRows) {
                return;
            }
            nRows = newNRows;
            if (rowHeaders != null) {
                rowHeaders = (String[])Arrays.resizeArray(rowHeaders, nRows);
            }
        }
        
        /**
         * Returns the number of rows.
         */
        public int getNRows() {
            return nRows;
        }
        
        /**
         * Returns a copy of the row headers or null if no row headers exist.
         */
        public String[] getRowHeaders() {
            return (rowHeaders == null) ? null : rowHeaders.clone();
        }
        
        /**
         * Sets the row headers.  null may be given to indicate no row headers.
         * If not null, nRows will be updated to match the number of row headers.  
         */
        public void setRowHeaders(String[] newRowHeaders) {
            rowHeaders = (newRowHeaders == null) ? null : newRowHeaders.clone();
            if (rowHeaders != null) {
                nRows = rowHeaders.length;
            }
        }
        
        /**
         * Returns a copy of the array of DataInfoFactories for the columns
         */
        public DataInfoDoubleArrayFactory[] getColumnInfoFactories() {
            return columnInfoFactories.clone();
        }
        
        /**
         * Sets the DataInfoFactories for the columns
         */
        public void setColumnInfoFactories(DataInfoDoubleArrayFactory[] newColumnInfoFactories) {
            columnInfoFactories = newColumnInfoFactories.clone();
        }
        
        private static final long serialVersionUID = 1L;
        protected DataInfoDoubleArrayFactory[] columnInfoFactories;
        protected String[] rowHeaders;
        protected int nRows;
    }
}
