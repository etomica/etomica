package etomica.data;

import etomica.data.types.DataDouble;
import etomica.data.types.DataGroup;
import etomica.data.types.DataGroup.DataInfoGroup;

/**
 * DataSink that takes a DataGroup object and pushes each sub-Data object from
 * the incoming DataGroup to a different stream.  The individual IDataInfo from
 * each subdata are pushed to each stream.  If a dataSink for one of the 
 * elements is not set or is set to null, the corresponding subdata element is
 * not pushed.
 *
 * @author Andrew Schultz
 */
public class DataGroupSplitter implements DataSink {

    public DataGroupSplitter() {
        tag = new DataTag();
        dataSinks = new DataSink[0];
    }
 
    /**
     * Returns the DataTag associated with this DataGroupSplitter
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

    public int getNumDataSinks() {
        return dataSinks.length;
    }
    
    /**
     * Sets the DataSink for the ith output stream (corresponding to the ith
     * numerical value coming in).
     */
    public void setDataSink(int i, DataSink newDataSink) {
        dataSinks[i] = newDataSink;
        if (dataSinks[i] != null && dataInfoGroup != null) {
            dataSinks[i].putDataInfo(dataInfoGroup.getSubDataInfo(i));
        }
    }

    public DataPipe getDataCaster(IDataInfo incomingDataInfo) {
        if (!(incomingDataInfo instanceof DataInfoGroup)) {
            throw new RuntimeException("I want to take a DataGroup");
        }
        return null;
    }

    public void putData(Data data) {
        DataGroup dataGroup = (DataGroup)data;
        for (int i=0; i<dataSinks.length; i++) {
            if (dataSinks[i] != null) {
                dataSinks[i].putData(dataGroup.getData(i));
            }
        }
    }

    public void putDataInfo(IDataInfo incomingDataInfo) {
        dataInfoGroup = (DataInfoGroup)incomingDataInfo;
        if (dataSinks.length != dataInfoGroup.getNDataInfo()) {
            dataSinks = new DataSink[dataInfoGroup.getNDataInfo()];
        }
        
        for (int i=0; i<dataSinks.length; i++) {
            if (dataSinks[i] != null) {
                dataSinks[i].putDataInfo(dataInfoGroup.getSubDataInfo(i));
            }
        }
    }

    protected DataSink[] dataSinks;
    protected DataDouble[] outData;
    protected DataInfoGroup dataInfoGroup;
    protected final DataTag tag;
}
