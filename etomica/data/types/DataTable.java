package etomica.data.types;

import java.io.Serializable;

import etomica.data.Data;
import etomica.data.DataFactory;
import etomica.data.DataInfo;
import etomica.units.Dimension;

/**
 * Data object that holds <tt>double[]</tt> arrays as if they are columns in a
 * table. The inner class <tt>Column</tt> wraps the array along with a
 * descriptive String heading and a Dimension instance that indicates the
 * physical dimensions of the data in the column.
 * <p>
 * The number of columns is set at construction and cannot be changed.
 * <p>
 * All columns are of equal length (i.e., the length of the <tt>double[]</tt>
 * arrays is the same for all columns).
 * 
 * @author David Kofke and Andrew Schultz
 *  
 */

/*
 * History Created on Jul 28, 2005 by kofke
 */
public class DataTable extends Data implements Serializable {

    /**
     * Creates a new table using the given arguments. The length of the DataInfo
     * array specifies the number of columns.
     * 
     * @param label
     *            a descriptive label for the table
     * @param columnDataInfo
     *            set of DataInfo instances used to set the heading and
     *            dimension of each column
     * @param columnLength
     *            the length of each and every column
     */
    public DataTable(String label, DataInfo[] columnDataInfo, int columnLength) {
        super(new DataInfo(label, evaluateDimension(columnDataInfo), new Factory(makeColumns(
                columnDataInfo, columnLength))));
        myColumns = ((DataTable.Factory) dataInfo.getDataFactory()).myColumns;
    }
    
    //determines if a given set of DataInfo are all of the same dimension;
    //if not the same, returns Dimension.MIXED
    //used by constructor above.
    private static Dimension evaluateDimension(DataInfo[] dataInfo) {
        Dimension dim = null;
        for(int i=0; i<dataInfo.length; i++) {
            if((dataInfo[i].getDimension() != dim) && (dim != null)) {
                return Dimension.MIXED;
            }
            dim = dataInfo[i].getDimension();
        }
        return dim;
    }

    /**
     * Creates a new table with a specified number of columns all of a given
     * length. The assigned heading for each column is its ordinal (starting
     * from 0), e.g., the heading of the first column is '0'.  All columns are
     * assigned the given Dimension, which is also the Dimension of the DataTable
     * (as kept in DataInfo).
     * 
     * @param label
     *            a descriptive label for the table
     * @param dimension
     *            indicates the common Dimension for all columns
     * @param numColumns
     *            the number of columns to make
     * @param columnLength
     *            the length of each and every column
     */
    public DataTable(String label, Dimension dimension, int numColumns,
            int columnLength) {
        super(new DataInfo(label, dimension, new Factory(makeColumns(
                dimension, numColumns, columnLength))));
        myColumns = ((DataTable.Factory) dataInfo.getDataFactory()).myColumns;
    }

    /**
     * Constructs a new table using the data in the given columns.
     * 
     * @param label
     *            a descriptive label for the table
     * @param myColumns
     *            the columns used to make the table (the instance is used
     *            direction, and is not copied/cloned)
     */
    public DataTable(String label, DataTable.Column[] myColumns) {
        super(new DataInfo(label, Dimension.UNDEFINED, new Factory(myColumns)));
        this.myColumns = myColumns;
    }

    //used by constructor
    private static Column[] makeColumns(DataInfo[] columnDataInfo,
            int columnLength) {
        Column[] columns = new Column[columnDataInfo.length];
        for (int i = 0; i < columnDataInfo.length; i++) {
            columns[i] = new Column(new double[columnLength], columnDataInfo[i]);
        }
        return columns;
    }

    //used by constructor
    private static Column[] makeColumns(Dimension dimension, int numColumns,
            int columnLength) {
        Column[] columns = new Column[numColumns];
        for (int i = 0; i < numColumns; i++) {
            columns[i] = new Column(new double[columnLength], Integer
                    .toString(i), dimension);
        }
        return columns;
    }

    /**
     * Copy constructor. Makes a new DataTable with new Column instances having
     * the same data values as in the given DataTable.
     */
    public DataTable(DataTable dataTable) {
        super(dataTable);
        myColumns = new DataTable.Column[dataTable.myColumns.length];
        for (int i = 0; i < myColumns.length; i++) {
            myColumns[i] = new DataTable.Column(dataTable.myColumns[i]);
        }
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
    public void E(Data data) {
        DataTable table = (DataTable) data;
        for (int i = 0; i < myColumns.length; i++) {
            System.arraycopy(table.myColumns[i].data, 0, myColumns[i].data, 0,
                    myColumns[i].data.length);
        }
    }

    /**
     * Returns the i-th column of the table, numbering from zero.
     */
    public Column getColumn(int i) {
        return myColumns[i];
    }

    /**
     * Returns the number of columns in the table.
     */
    public int getNColumns() {
        return myColumns.length;
    }

    /**
     * Returns the number of rows in each and every column of the table.
     */
    public int getNRows() {
        return myColumns[0].data.length;
    }

    final DataTable.Column[] myColumns;

    /**
     * The inner-class Factory is configured to make another DataTable of the
     * same shape (nCols, nRows) as this one.
     */
    public static class Factory implements DataFactory, Serializable {

        final DataTable.Column[] myColumns;
        String[] rowHeaders;
        String[] columnHeaders;
        int columnLength;

        Factory(DataTable.Column[] myColumns) {
            this.myColumns = myColumns;
        }

        /**
         * Makes a new DataTable, assigning the descriptive label to it; the
         * given Dimension is not used, but is included as an argument to adhere
         * to the DataFactory interface.
         */
        public Data makeData(String label, Dimension dimension) {
            DataTable.Column[] newColumns = new DataTable.Column[myColumns.length];
            for (int i = 0; i < newColumns.length; i++) {
                newColumns[i] = new DataTable.Column(myColumns[i]);
            }
            //drops dimension
            return new DataTable(label, newColumns);
        }

        /**
         * Returnsw DataTable.class.
         */
        public Class getDataClass() {
            return DataTable.class;
        }

        /**
         * Returns the number of columns in the DataTable made by this
         * DataFactory.
         * 
         * @return
         */
        public int getNColumns() {
            return columnHeaders.length;
        }

        /**
         * Returns the length of each and every row in the DataTable made by
         * this DataFactory.
         */
        public int getNRows() {
            return columnLength;
        }
    }

    /**
     * Class that wraps a <tt>double[]</tt> array, a descriptive String that
     * is suitable as a heading for the column, and a Dimension indicating the
     * physical dimensions of the data. The wrapped data array is final, so the
     * length of the column cannot be changed after construction.
     */
    public static class Column implements Serializable {

        /**
         * Makes a Column wrapping the given array and using the DataInfo to
         * determine the heading and dimension.
         */
        public Column(double[] data, DataInfo dataInfo) {
            this(data, dataInfo.getLabel(), dataInfo.getDimension());
        }

        /**
         * Makes a Column wrapping the given array and using the given values
         * for the heading and dimension.
         */
        public Column(double[] data, String label, Dimension dimension) {
            this.data = data;
            this.heading = label;
            this.dimension = dimension;
        }

        /**
         * Copy constructor. Makes independent copies of all data in column.
         */
        public Column(Column column) {
            this.data = (double[]) column.data.clone();
            this.heading = column.heading;
            this.dimension = column.dimension;
        }

        /**
         * Mutator method for the heading.
         */
        public void setHeading(String heading) {
            this.heading = heading;
        }

        /**
         * Accessor method for the heading.
         */
        public String getHeading() {
            return heading;
        }

        /**
         * Accessor method for the (immutable) dimension.
         */
        public Dimension getDimension() {
            return dimension;
        }

        /**
         * Returns the wrapped array.
         */
        public double[] getData() {
            return data;
        }

        private final double[] data;
        private String heading;
        private final Dimension dimension;
    }

}