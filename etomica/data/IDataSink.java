package etomica.data;

import etomica.api.IData;


/**
 * A recipient of Data.  Data goes in and might (or might not) come out.
 */
public interface IDataSink {

    /**
     * Gives data to DataSink for processing, display, or whatever it does.
     */
    public void putData(IData data);

    /**
     * Informs the DataSink of the type of Data it should expect to receive.
     */
    public void putDataInfo(IEtomicaDataInfo dataInfo);

    /**
     * Returns a DataProcessor that casts the data that will be given 
     * to this DataSink to a form that it can accept.  Returns null if the
     * Data is acceptable without casting.  Otherwise object that is pushing
     * the Data into this DataSink is responsible for filtering the Data through
     * the returned DataTransformer. 
     * 
     * @param dataInfo the DataInfo for the Data that will fed to the sink's putData method
     */
    public DataPipe getDataCaster(IEtomicaDataInfo dataInfo);
}