package etomica.data;

import etomica.Action;
import etomica.DataPipe;
import etomica.DataSink;
import etomica.DataSource;

/**
 * A data pipe that can move data from a source to its sinks
 * as the result of an action performed.
 */
public class DataPump extends DataPipe implements Action {

	/**
	 * Constructs DataPump with the given DataSource and
	 * DataSinks.  Data source cannot be null.
	 */
	public DataPump(DataSource dataSource, DataSink[] dataSinks) {
		if(dataSource == null) throw new NullPointerException("Error: cannot construct data pump without a data source");
		this.dataSource = dataSource;
		setDataSinks(dataSinks);
	}
    
    /**
     * Constructs DataAccumulator with the given DataSource and a single DataSink.
     */
    public DataPump(DataSource dataSource, DataSink dataSink) {
        this(dataSource, (dataSink == null) ? null : new DataSink[] {dataSink});
    }
	
	/**
     * Transmits the data from the source to all sinks.
	 */
	public void actionPerformed() {
        pushData(dataSource.getData());
    }
    
    /**
     * Passes the given data on to all sinks.  Normally invoked
     * via actionPerformed, which gets the data from the data source
     * specified at construction.  If this method is called directly,
     * the data pump behaves as a simple data pipe.
     */
    public void putData(double[] data) {
        pushData(data);
	}
    	
    /**
     * @return Returns the dataSource.
     */
    public DataSource getDataSource() {
        return dataSource;
    }

    private final DataSource dataSource;
    
}
