package etomica.data;

import etomica.Data;
import etomica.DataInfo;
import etomica.DataSink;

/**
 * An object that receives Data, processes it, and pushes the result on to a DataSink.
 */

public abstract class DataProcessor implements DataPipe, java.io.Serializable {

    /**
     * Processes the input Data and returns Data for pushing to the next DataSink.
     * Returns null if output Data should not (yet) be pushed downstream.
     * @param inputData the Data for processing
     * @return the processed Data, for sending downstream (if not null)
     */
    protected abstract Data processData(Data inputData);
    
    /**
     * Returns the DataInfo of the Data that is output from this DataProcessor.  This
     * also notifies the subclass of this information, so that it can make any
     * objects or otherwise prepare processData.
     * 
     * @inputDataInfo the DataInfo of the Data that will be input to this DataProcessor
     * @return the DataInfo of the Data that will be output by this DataProcessor
     */
    protected abstract DataInfo processDataInfo(DataInfo inputDataInfo);
    

    /**
     * Processes input Data and pushes it downstream if output Data and
     * DataSink are not null. 
     */
    public void putData(Data data) {
        Data outputData = processData(data);
        if(dataSink != null && outputData != null) {
            dataSink.putData(outputData);
        }
    }
    
    public void putDataInfo(DataInfo inputDataInfo) {
        dataInfo = processDataInfo(inputDataInfo);
        insertTransformerIfNeeded();
        if(dataSink != null) {
            dataSink.putDataInfo(dataInfo);
        }
    }

    /**
     * @return Returns the data sink.
     */
    public DataSink getDataSink() {
        return dataSink;
    }

    /**
     * Sets the sink receiving the data. Null value is permitted.
     * 
     * @param dataSinks
     *            The data sinks to set.
     */
    public void setDataSink(DataSink dataSink) {
        this.dataSink = dataSink;
        insertedTransformer = false;
        insertTransformerIfNeeded();
        if(dataSink != null && dataInfo != null) {
            dataSink.putDataInfo(dataInfo);
        }
    }
    
    private void insertTransformerIfNeeded() {
        if(dataSink == null || dataInfo == null) return;
        //remove transformer if one was previously inserted
        if(insertedTransformer) {
            dataSink = ((DataProcessor)dataSink).dataSink;
            insertedTransformer = false;
        }
        DataProcessor caster = dataSink.getDataCaster(dataInfo);
        if(caster != null) {
            caster.setDataSink(dataSink);
            dataSink = caster;
            insertedTransformer = true;
        }
    }

    protected DataSink dataSink;
    protected DataInfo dataInfo;
    private boolean insertedTransformer = false;

}
