/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.data;


/**
 * An object that receives Data, processes it, and pushes the result on to a
 * DataSink.
 */

public abstract class DataProcessor implements DataPipe {

    public DataProcessor() {
        tag = new DataTag();
    }
    
    public DataTag getTag() {
        return tag;
    }
    
    /**
     * Processes the input Data and returns Data for pushing to the next
     * DataSink. Returns null if output Data should not (yet) be pushed
     * downstream.
     * 
     * @param inputData
     *            the Data for processing
     * @return the processed Data, for sending downstream (if not null)
     */
    protected abstract IData processData(IData inputData);

    /**
     * Informs this DataProcessor of the DataInfo for the Data it will be
     * processing. Typically the subclass will use this information to make any
     * objects or otherwise prepare for calls to processData.
     * 
     * @param inputDataInfo the DataInfo of the Data that will be input to this
     *                DataProcessor
     * @return the DataInfo of the Data that will be output by this
     *         DataProcessor
     */
    protected abstract IEtomicaDataInfo processDataInfo(IEtomicaDataInfo inputDataInfo);

    /**
     * Processes input Data and pushes it downstream if output Data and DataSink
     * are not null.
     */
    public void putData(IData data) {
        IData outputData = processData(data);
        if (dataSink != null && outputData != null) {
            dataSink.putData(outputData);
        }
    }

    /**
     * Invokes processDataInfo on the given DataInfo, and passes the returned
     * DataInfo to the dataSink (if not null).  Will insert a data caster before
     * the DataSink if appropriate.
     */
    public void putDataInfo(IEtomicaDataInfo inputDataInfo) {
        dataInfo = processDataInfo(inputDataInfo);
        insertTransformerIfNeeded();
        if (dataSink != null && dataInfo != null) {
            dataSink.putDataInfo(dataInfo);
        }
    }
    
    public IEtomicaDataInfo getDataInfo() {
        return dataInfo;
    }

    /**
     * @return Returns the data sink, which may be null.
     */
    public IDataSink getDataSink() {
        return trueDataSink;
    }

    /**
     * Sets the sink receiving the data. Null value is permitted.
     * 
     * @param newDataSink
     *            The data sink to set.
     */
    public void setDataSink(IDataSink newDataSink) {
        //trueDataSink is the sink that the caller acutally cares about
        //dataSink is the immeadiate sink for this processor (might be a transformer)
        trueDataSink = newDataSink;
        insertTransformerIfNeeded();
        if (dataSink != null && dataInfo != null) {
            dataSink.putDataInfo(dataInfo);
        }
    }

    protected void insertTransformerIfNeeded() {
        if (trueDataSink == null || dataInfo == null)
            return;
        //remove transformer if one was previously inserted
        dataSink = trueDataSink;
        DataPipe caster = dataSink.getDataCaster(dataInfo);
        if (caster != null) {
            caster.setDataSink(trueDataSink);
            dataSink = caster;
        }
    }

    protected IDataSink dataSink;
    protected IDataSink trueDataSink;
    protected IEtomicaDataInfo dataInfo;
    protected final DataTag tag;
}
