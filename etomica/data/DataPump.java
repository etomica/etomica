package etomica.data;

import etomica.Action;
import etomica.DataSink;
import etomica.DataSource;
import etomica.DataTranslator;

/**
 * A data pusher that can move data from a source to its sinks
 * as the result of an action performed.
 */
public class DataPump extends DataPusher implements Action {

	/**
	 * Constructs DataPump with the given DataSource and
	 * DataSinks.  Data source cannot be null.
	 */
	public DataPump(DataSource dataSource, DataSink[] dataSinks) {
		if(dataSource == null) throw new NullPointerException("Error: cannot construct data pump without a data source");
		this.dataSource = dataSource;
        setDefaultLabel(dataSource.getLabel());
        setDimension(dataSource.getDimension());
		setDataSinks(dataSinks);
	}
    
    /**
     * Constructs DataPump with the given DataSource and a single DataSink.
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
     * @return Returns the dataSource.
     */
    public DataSource getDataSource() {
        return dataSource;
    }

    public DataTranslator getTranslator() {
        return DataTranslator.IDENTITY;
    }
    
    private final DataSource dataSource;
    
}
