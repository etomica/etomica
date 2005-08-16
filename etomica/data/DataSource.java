package etomica.data;

/**
 * Interface for an object that can provide Data objects
 * on request.  Normally a DataSource heads a stream that processes
 * and/or records the Data as it passes through different segments.  
 * Data is pulled from the DataSource by a DataPump and pushed down the
 * data stream.
 */
 
public interface DataSource {

    /**
     * @return the data given by this source
     */
	public Data getData();

    /**
     * Returns the DataInfo instance that will be held by Data
     * given by this source.  This information is useful for
     * setting up the data stream and for providing annotation when
     * displaying or writing the Data.
     */
    public DataInfo getDataInfo();
    
}