/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

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
public class DataGroupSplitter implements IDataSink {

    protected final DataTag tag;
    protected IDataSink[] dataSinks;
    protected DataDouble[] outData;
    protected DataInfoGroup dataInfoGroup;
    
    public DataGroupSplitter() {
        tag = new DataTag();
        dataSinks = new IDataSink[0];
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
    public IDataSink getDataSink(int i) {
        return dataSinks[i];
    }

    public int getNumDataSinks() {
        return dataSinks.length;
    }

    /**
     * Sets the DataSink for the ith output stream (corresponding to the ith
     * numerical value coming in).
     */
    public void setDataSink(int i, IDataSink newDataSink) {
        dataSinks[i] = newDataSink;
        if (dataSinks[i] != null && dataInfoGroup != null) {
            dataSinks[i].putDataInfo(dataInfoGroup.getSubDataInfo(i));
        }
    }

    public void putData(IData data) {
        DataGroup dataGroup = (DataGroup)data;
        for (int i=0; i<dataSinks.length; i++) {
            if (dataSinks[i] != null) {
                dataSinks[i].putData(dataGroup.getData(i));
            }
        }
    }

    public void putDataInfo(IDataInfo incomingDataInfo) {
        if (!(incomingDataInfo instanceof DataInfoGroup)) {
            throw new RuntimeException("I want to take a DataGroup");
        }
        dataInfoGroup = (DataInfoGroup)incomingDataInfo;
        if (dataSinks.length != dataInfoGroup.getNDataInfo()) {
            dataSinks = new IDataSink[dataInfoGroup.getNDataInfo()];
        }

        for (int i=0; i<dataSinks.length; i++) {
            if (dataSinks[i] != null) {
                dataSinks[i].putDataInfo(dataInfoGroup.getSubDataInfo(i));
            }
        }
    }
}
