package etomica.data;

/**
 * Interface for an object through which data can flow. It can serve as a
 * DataSink and pass Data onto another DataSink.
 * 
 * @author David Kofke
 *  
 */
public interface DataPipe extends IDataSink {

    /**
     * Sets the DataSink that will receive data from this DataPipe. Any
     * previously set DataSink(s) is discarded.
     * 
     * @param dataSink
     *            the new recipient of Data from this DataPipe.
     */
    public void setDataSink(IDataSink dataSink);
    
    public DataTag getTag();

}
