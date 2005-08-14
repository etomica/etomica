package etomica.data.types;

import java.io.Serializable;

import etomica.data.Data;
import etomica.data.DataFactory;
import etomica.data.DataInfo;
import etomica.units.Dimension;



/*
 * History
 * Created on Jul 28, 2005 by kofke
 */
/**
 * Data object that is sent from the data caster to the table.  Holds the columns
 * that were updated.
 */
public class DataTable extends Data implements Serializable {

    public DataTable(String label, DataInfo[] columnDataInfo, int columnLength) {
        super(new DataInfo(label, Dimension.UNDEFINED, new Factory(makeColumns(columnDataInfo, columnLength))));
        myColumns = ((DataTable.Factory)dataInfo.getDataFactory()).myColumns;
        index = 5;
    }
    
    public DataTable(String label, Dimension dimension, int numColumns, int columnLength) {
        super(new DataInfo(label, Dimension.UNDEFINED, new Factory(makeColumns(dimension,numColumns,columnLength))));
        myColumns = ((DataTable.Factory)dataInfo.getDataFactory()).myColumns;
        index = 5;
    }

    public DataTable(int index, String label, DataTable.Column[] myColumns) {
        super(new DataInfo(label, Dimension.UNDEFINED, new Factory(myColumns)));
        this.index = index;
        this.myColumns = myColumns;
    }
    
    private static Column[] makeColumns(DataInfo[] columnDataInfo, int columnLength) {
        Column[] columns = new Column[columnDataInfo.length];
        for (int i=0; i<columnDataInfo.length; i++) {
            columns[i] = new Column(new double[columnLength],columnDataInfo[i]);
        }
        return columns;
    }
    
    private static Column[] makeColumns(Dimension dimension, int numColumns, int columnLength) {
        Column[] columns = new Column[numColumns];
        for (int i=0; i<numColumns; i++) {
            columns[i] = new Column(new double[columnLength],Integer.toString(i),dimension);
        }
        return columns;
    }
    
    public DataTable(DataTable dataTable) {
        super(dataTable);
        this.index = dataTable.index;
        myColumns = new DataTable.Column[dataTable.myColumns.length];
        for(int i=0; i<myColumns.length; i++) {
            myColumns[i] = new DataTable.Column(dataTable.myColumns[i]);
        }
    }
    
    public Data makeCopy() {
        return new DataTable(this);
    }
    
    //copies only data, not headings
    public void E(Data data) {
        DataTable table = (DataTable)data;
        for(int i=0; i<myColumns.length; i++) {
            System.arraycopy(table.myColumns[i].data, 0, myColumns[i].data, 0, myColumns[i].data.length);
        }
    }
    
    public Column getColumn(int i) {
        return myColumns[i];
    }
    
    public int getNColumns() {
        return myColumns.length;
    }
    
    public int getNRows() {
        return myColumns[0].data.length;
    }

    final int index;
    final DataTable.Column[] myColumns;
    
    public static class Factory implements DataFactory, Serializable {

        final DataTable.Column[] myColumns;
        String[] rowHeaders;
        String[] columnHeaders;
        int columnLength;
        
        Factory(DataTable.Column[] myColumns) {
            this.myColumns = myColumns;
        }
        public Data makeData(String label, Dimension dimension) {
            DataTable.Column[] newColumns = new DataTable.Column[myColumns.length];
            for(int i=0; i<newColumns.length; i++) {
                newColumns[i] = new DataTable.Column(myColumns[i]);
            }
            //drops dimension
            return new DataTable(-1, label, newColumns);
        }
        
        public Class getDataClass() {
            return DataTable.class;
        }
        
        public int getNColumns() {
            return columnHeaders.length;
        }
        
        public int getNRows() {
            return columnLength;
        }
    }
    
    public static class Column implements Serializable {
        
        public Column(double[] data, DataInfo dataInfo) {
            this(data,dataInfo.getLabel(),dataInfo.getDimension());
        }
        public Column(double[] data, String label, Dimension dimension) {
            this.data = data;
            this.heading = label;
            this.dimension = dimension;
        }
        /**
         * Copy constructor.  Makes independent copies of all data in column.
         */
        public Column(Column column) {
            this.data = (double[])column.data.clone();
            this.heading = column.heading;
            this.dimension = column.dimension;
        }
        
        public void setHeading(String heading) {
            this.heading = heading;
        }
        
        public String getHeading() {
            return heading;
        }
        
        public Dimension getDimension() {
            return dimension;
        }
        
        public double[] getData() {
            return data;
        }
        
        private final double[] data;
        private String heading;
        private final Dimension dimension;
    }

}