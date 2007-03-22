package etomica.data;

import etomica.data.types.DataDouble;

/**
 * DataSink that takes a Data object and pushes each numerical value from the
 * incoming Data to a different stream.  This is not a DataPipe because
 * setDataSink would be inappropriate.
 *
 * @author Andrew Schultz
 */
public class DataSplitter implements DataSink {

    public DataSplitter() {
        tag = new DataTag();
        dataSinks = new DataSink[0];
    }
 
    /**
     * Returns the DataTag associated with this DataSplitter
     * @return
     */
    public DataTag getTag() {
        return tag;
    }

    /**
     * Returns the DataSink for the ith output stream (corresponding to the ith
     * numerical value coming in).
     */
    public DataSink getDataSink(int i) {
        return dataSinks[i];
    }

    /**
     * Sets the DataSink for the ith output stream (corresponding to the ith
     * numerical value coming in).
     */
    public void setDataSink(int i, DataSink newDataSink) {
        dataSinks[i] = newDataSink;
        if (dataSinks[i] != null && dataInfo != null) {
            dataSinks[i].putDataInfo(dataInfo);
        }
    }

    public DataPipe getDataCaster(IDataInfo incomingDataInfo) {
        return null;
    }

    public void putData(Data data) {
        for (int i=0; i<dataSinks.length; i++) {
            if (dataSinks[i] != null) {
                outData[i].x = data.getValue(i);
                dataSinks[i].putData(outData[i]);
            }
        }
    }

    public void putDataInfo(IDataInfo incomingDataInfo) {
        if (dataSinks.length != incomingDataInfo.getLength()) {
            dataSinks = new DataSink[incomingDataInfo.getLength()];
            
            //do we really need a separate out data for each value?
            outData = new DataDouble[incomingDataInfo.getLength()];
            for (int i=0; i<outData.length; i++) {
                outData[i] = new DataDouble();
            }
        }
        
        dataInfo = new DataDouble.DataInfoDouble(incomingDataInfo.getLabel(), incomingDataInfo.getDimension());
        
        for (int i=0; i<dataSinks.length; i++) {
            if (dataSinks[i] != null) {
                dataSinks[i].putDataInfo(dataInfo);
            }
        }
    }

    protected DataSink[] dataSinks;
    protected DataDouble[] outData;
    protected DataInfo dataInfo;
    protected final DataTag tag;
}
