package etomica.data;

import etomica.Action;
import etomica.Data;
import etomica.DataInfo;
import etomica.DataSink;
import etomica.DataSource;

/**
 * A data pusher that can move data from a source to its sinks
 * as the result of an action performed.
 */
public class DataPump extends DataProcessor implements Action {

	/**
	 * Constructs DataPump with the given DataSource and
	 * DataSink.  Data source cannot be null.
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
     * Transmits the data from the source to all sinks.
	 */
	public void actionPerformed() {
        Data data = dataSource.getData();
        if (dataSourceInfo != data.getDataInfo()) {
            dataSourceInfo = data.getDataInfo();
            if (dataSink != null) {
                dataSink.putDataInfo(dataSourceInfo);
            }
        }
        putData(data);
    }
    
    public Data processData(Data inputData) {
        return inputData;
    }
    
    public DataInfo processDataInfo(DataInfo inputDataInfo) {
        return inputDataInfo;
    }
    
    /**
     * Returns null, indicating that this DataSink can handle any type of Data without casting.
     */
    public DataProcessor getDataCaster(DataInfo dataInfo) {
        return null;
    }
    
    /**
     * @return Returns the dataSource.
     */
    public DataSource getDataSource() {
        return dataSource;
    }

    public String getLabel() {
        return label;
    }
    
    public void setLabel(String label) {
        this.label = label;
    }
    
    private DataInfo dataSourceInfo;
    private final DataSource dataSource;
    private String label;
}
