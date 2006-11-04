package etomica.data;

import etomica.action.Action;

/**
 * A DataProcessor whose action is to actively take Data from a DataSource and send it to
 * DataSinks.  
 */
public class DataPump extends DataProcessor implements Action {

	/**
	 * Constructs DataPump with the given DataSource and
	 * DataSink.  Data source cannot be null.  Data sink can 
     * be null and must be identified via setDataSink if DataPump
     * is to have any effect.
	 */
    public DataPump(DataSource dataSource, DataSink dataSink) {
        if(dataSource == null) throw new NullPointerException("Error: cannot construct data pump without a data source");
        this.dataSource = dataSource;
        dataSourceInfo = dataSource.getDataInfo();
        setDataSink(dataSink);
        setLabel("Data Pump");
        putDataInfo(dataSource.getDataInfo());
	}
    
	/**
     * Transmits the data from the source to the sink. Before transmitting
     * the Data, this method will first check that the DataInfo from the source
     * is the same as it was last time this method was invoked.  If it has changed,
     * a call to putDataInfo in the sink will be invoked before passing along the Data.
	 */
	public void actionPerformed() {
        Data data = dataSource.getData();
        if (dataSourceInfo != dataSource.getDataInfo()) {
            dataSourceInfo = dataSource.getDataInfo();
            if (dataSink != null) {
                dataSink.putDataInfo(dataSourceInfo);
            }
        }
        putData(data);
    }
    
    /**
     * Returns the given Data.
     */
    public Data processData(Data inputData) {
        return inputData;
    }
    
    /**
     * Returns the given DataInfo.
     */
    public DataInfo processDataInfo(DataInfo inputDataInfo) {
        dataInfo = inputDataInfo.getFactory().makeDataInfo();
        dataInfo.addTag(getTag());
        return dataInfo;
    }
    
    /**
     * Returns null, indicating that this DataSink can handle any type of Data without casting.
     */
    public DataProcessor getDataCaster(DataInfo incomingDataInfo) {
        return null;
    }
    
    /**
     * @return Returns the dataSource.
     */
    public DataSource getDataSource() {
        return dataSource;
    }

    /**
     * Returns a descriptive label for this DataPump.  Part of the Action interface.
     * Has no connection to the label for the Data going through.
     */
    public String getLabel() {
        return label;
    }
    
    /**
     * Sets a descriptive label for this DataPump.
     */
    public void setLabel(String label) {
        this.label = label;
    }
    
    private DataInfo dataSourceInfo;
    private final DataSource dataSource;
    private String label;
}
