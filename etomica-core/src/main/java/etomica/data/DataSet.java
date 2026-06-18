/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.data;

import etomica.data.types.DataGroup;
import etomica.data.types.DataGroup.DataInfoGroup;
import etomica.graphics.DisplayPlot;
import etomica.graphics.DisplayPlotXChart;
import etomica.util.Arrays;

/**
 * Receives and organizes data from multiple streams. Fires event when new
 * data arrives; can be configured to fire when any piece of data
 * changes, or to perform action only when all data has changed.
 * 
 * @author David Kofke
 */

public class DataSet {
    private static int SINK_COUNT = 0;

    protected DataSetListener[] listeners = new DataSetListener[0];
    protected DataSetSink[] psuedoSinks;
    protected int[] forwardDataMap;
    protected int[] backwardDataMap;
    private boolean updatingOnAnyChange;
    private int[] dataChangedBits = new int[0];

    public DataSet() {
        psuedoSinks = new DataSetSink[0];
        backwardDataMap = new int[0];
        forwardDataMap = new int[0];
    }

    public DataSetSink makeDataSink() {
        if (this.listeners.length > 0 && this.listeners[0] instanceof DisplayPlotXChart) {
            return (DataSetSink) ((DisplayPlotXChart) this.listeners[0]).makeSink("Sink" + SINK_COUNT++);
        }
        return new DataSetSink(this);
    }

    public DataSetSink makeDataSink(String name) {
        return new DataSetSink(this, name);
    }

    public void reset() {
        for (int i=0; i<psuedoSinks.length; i++) {
            psuedoSinks[i].index = -1;
        }
        psuedoSinks = new DataSetSink[0];
        backwardDataMap = new int[0];
        forwardDataMap = new int[0];
    }

    /**
     * Returns the ith Data from the set.
     */
    public IData getData(int i) {
        int iData = backwardDataMap[i];
        IData data = psuedoSinks[iData].getData();
        if (data instanceof DataGroup) {
            data = ((DataGroup)data).getData(i-forwardDataMap[iData]);
        }
        return data;
    }

    /**
     * Returns the ith DataInfo from the set.
     */
    public IDataInfo getDataInfo(int i) {
        int iData = backwardDataMap[i];
        IDataInfo dataInfo = psuedoSinks[iData].getDataInfo();
        if (dataInfo instanceof DataInfoGroup) {
            dataInfo = ((DataInfoGroup)dataInfo).getSubDataInfo(i-forwardDataMap[iData]);
        }
        return dataInfo;
    }

    public String getName(int i) {
        int iData = backwardDataMap[i];
        return psuedoSinks[iData].getName();
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
                dataChangedBits = java.util.Arrays.copyOf(dataChangedBits,dataChangedBits.length+1);
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
            IDataInfo dataInfo = psuedoSinks[i].getDataInfo();
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
    
    /**
     * A special DataSink used by the DataSet to receive Data and DataInfo 
     * from a single stream.
     * 
     * @author Andrew Schultz
     */
    protected static class DataSetSink implements IDataSink {
        private final DataSet dataSet;
        // index exists solely to indicate whether this DataSetSink is not in
        // the array of DataSetSinks (index=-1) or its position in the array so
        // the DataSet doesn't have to go looking to find it.
        protected int index;
        private IDataInfo incomingDataInfo;
        private IData data;
        private final String name;

        public DataSetSink(DataSet dataSet, String name) {
            this.dataSet = dataSet;
            this.name = name;
            this.index = -1;
        }

        public DataSetSink(DataSet dataSet) {
            this(dataSet, "Sink" + SINK_COUNT++);
        }

        public String getName() {
            return this.name;
        }
        
        public void putDataInfo(IDataInfo newDataInfo) {
            incomingDataInfo = newDataInfo;
            dataSet.dataInfoChanged(this);
        }

        public void putData(IData incomingData) {
            data = incomingData;
            dataSet.dataChanged(this);
        }

        public IDataInfo getDataInfo() {
            return incomingDataInfo;
        }

        public IData getData() {
            return data;
        }
    }
}
