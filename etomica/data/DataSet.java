package etomica.data;

import java.io.Serializable;

import etomica.data.types.DataGroup;
import etomica.data.types.DataGroup.DataInfoGroup;
import etomica.util.Arrays;

/**
 * Receives and organizes data from multiple streams. Fires event when new
 * data arrives; can be configured to fire when any piece of data
 * changes, or to perform action only when all data has changed.
 * 
 * @author David Kofke
 */

public class DataSet implements Serializable {

    public DataSet() {
        this(null);
    }
    
    public DataSet(DataCasterJudge dataCasterJudge) {
        psuedoSinks = new DataSetSink[0];
        this.dataCasterJudge = dataCasterJudge;
        backwardDataMap = new int[0];
        forwardDataMap = new int[0];
    }
    
    public DataCasterJudge getDataCasterJudge() {
        return dataCasterJudge;
    }

    public DataSetSink makeDataSink() {
        return new DataSetSink(dataCasterJudge, this);
    }
    
    public void reset() {
        for (int i=0; i<psuedoSinks.length; i++) {
            psuedoSinks[i].index = -1;
        }
        psuedoSinks = new DataSetSink[0];
    }
    
    /**
     * Returns the ith Data from the set.
     */
    public Data getData(int i) {
        int iData = backwardDataMap[i];
        Data data = psuedoSinks[iData].getData();
        if (data instanceof DataGroup) {
            data = ((DataGroup)data).getData(i-forwardDataMap[iData]);
        }
        return data;
    }        
    
    /**
     * Returns the ith DataInfo from the set.
     */
    public DataInfo getDataInfo(int i) {
        int iData = backwardDataMap[i];
        DataInfo dataInfo = psuedoSinks[iData].getDataInfo();
        if (dataInfo instanceof DataInfoGroup) {
            dataInfo = ((DataInfoGroup)dataInfo).getSubDataInfo(i-forwardDataMap[iData]);
        }
        return dataInfo;
    }

    public int getDataCount() {
        return backwardDataMap.length;
    }

    /**
     * @return flag defining protocol for firing table events that
     * notify of column data changes.
     */
    public boolean isUpdatingOnAnyChange() {
        return updatingOnAnyChange;
    }

    /**
     * Describes protocol for firing table events that notify of column data
     * changes. If true, event is fired when any column notifies that it has
     * changed; if false, event is fired only when all columns have
     * notified that they have changed (since last time event was fired).
     */
    public void setUpdatingOnAnyChange(boolean updatingWithAnyChange) {
        this.updatingOnAnyChange = updatingWithAnyChange;
    }

    /**
     * Notifies the DataSet that new Data has arrived for the given psuedo 
     * DataSink.
     */
    protected void dataChanged(DataSetSink dataSetSink) {
        boolean newData = false;
        int index = dataSetSink.index;
        if (index == -1) {
            newData = true;
            index = psuedoSinks.length; 

            dataSetSink.index = index;
            psuedoSinks = (DataSetSink[])Arrays.addObject(psuedoSinks, dataSetSink);

            if (index/32 > dataChangedBits.length-1) {
                dataChangedBits = Arrays.resizeArray(dataChangedBits,dataChangedBits.length+1);
            }
            dataChangedBits[dataChangedBits.length-1] |= 1<<(index%32);

            updateDataMap();
            fireDataCountChangedEvent();
        }

        // unset the bit to indicate this data object has been updated
        dataChangedBits[index >> 5] &= (0xFFFFFFFF ^ (1<<(index & 0xFFFFFFFF)));

        if (updatingOnAnyChange) {
            fireDataChangedEvent();
        }
        else {
            int l = dataChangedBits.length;
            // check to see if all data has been updated
            // always redraw when new data comes in
            if (!newData) {
                for (int i=0; i<l; i++) {
                    if (dataChangedBits[i] != 0) {
                        return;
                    }
                }
            }
            fireDataChangedEvent();
            // reset all bits but the last one
            for (int i=0; i<l-1; i++) {
                dataChangedBits[i] = 0xFFFFFFFF;
            }
            if (index%32 == 31) {
                // reset the last one
                dataChangedBits[l-1] = 0xFFFFFFFF;
            }
            else {
                // partially reset the last one
                dataChangedBits[l-1] = (1<<(index%32));
            }
        }
    }

    /**
     * Notifies the DataSet that new DataInfo has arrived for the given 
     * psuedo DataSink.
     */
    protected void dataInfoChanged(DataSetSink dataSetSink) {
        if (dataSetSink.index != -1) {
            // we've already seen this data, so it's in the map.
            updateDataMap();
        }
    }

    /**
     * This updates the mapping between pseudo DataSinks and indices (used for
     * getData(int), etc.).
     */
    protected void updateDataMap() {
        //forward maps the ith pseudo DataSink to the appropriate Data element
        forwardDataMap = new int[psuedoSinks.length];
        int dataCount = 0;
        for (int i=0; i<psuedoSinks.length; i++) {
            forwardDataMap[i] = dataCount;
            DataInfo dataInfo = psuedoSinks[i].getDataInfo();
            if (dataInfo instanceof DataInfoGroup) {
                dataCount += ((DataInfoGroup)dataInfo).getNDataInfo();
            }
            else {
                dataCount++;
            }
        }
        
        //backward maps the ith Data element to the appropriate psuedoDataSink
        backwardDataMap = new int[dataCount];
        for (int i=0; i<psuedoSinks.length-1; i++) {
            for (int j=forwardDataMap[i]; j<forwardDataMap[i+1]; j++) {
                backwardDataMap[j] = i;
            }
        }
        for (int j=forwardDataMap[psuedoSinks.length-1]; j<dataCount; j++) {
            backwardDataMap[j] = forwardDataMap.length-1;
        }
    }

    protected void fireDataChangedEvent() {
        for (int i = 0; i < listeners.length; i++) {
            listeners[i].dataChanged(this);
        }
    }

    protected void fireDataCountChangedEvent() {
        for (int i = 0; i < listeners.length; i++) {
            listeners[i].dataCountChanged(this);
        }
    }

    public void addDataListener(DataSetListener newListener) {
        listeners = (DataSetListener[]) Arrays.addObject(listeners,
                newListener);
    }

    public void removeDataListener(DataSetListener oldListener) {
        listeners = (DataSetListener[]) Arrays.removeObject(listeners,
                oldListener);
    }

    private static final long serialVersionUID = 1L;
    protected DataSetListener[] listeners = new DataSetListener[0];
    private boolean updatingOnAnyChange;
    private int[] dataChangedBits = new int[0];
    private final DataCasterJudge dataCasterJudge;
    protected DataSetSink[] psuedoSinks;
    
    protected int[] forwardDataMap;
    protected int[] backwardDataMap;
    
    public interface DataCasterJudge {
        public DataProcessor getDataCaster(DataInfo inputDataInfo);
    }
    
    /**
     * A special DataSink used by the DataSet to receive Data and DataInfo 
     * from a single stream.
     * 
     * @author Andrew Schultz
     */
    protected static class DataSetSink implements DataSink {
        public DataSetSink(DataCasterJudge dataCasterJudge, DataSet dataSet) {
            this.dataCasterJudge = dataCasterJudge;
            this.dataSet = dataSet;
            index = -1;
        }
        
        public void putDataInfo(DataInfo newDataInfo) {
            incomingDataInfo = newDataInfo;
            dataSet.dataInfoChanged(this);
        }
        
        public void putData(Data incomingData) {
            data = incomingData;
            dataSet.dataChanged(this);
        }
        
        public DataProcessor getDataCaster(DataInfo newIncomingDataInfo) {
            if (dataCasterJudge == null) {
                return null;
            }
            return dataCasterJudge.getDataCaster(newIncomingDataInfo);
        }
        
        public DataInfo getDataInfo() {
            return incomingDataInfo;
        }
        
        public Data getData() {
            return data;
        }

        private static final long serialVersionUID = 1L;
        private final DataCasterJudge dataCasterJudge;
        private DataInfo incomingDataInfo;
        private final DataSet dataSet;
        private Data data;
        // index exists solely to indicate whether this DataSetSink is not in
        // the array of DataSetSinks (index=-1) or its position in the array so
        // the DataSet doesn't have to go looking to find it.
        protected int index;
    }
}
