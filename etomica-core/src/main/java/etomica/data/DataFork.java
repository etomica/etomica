/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.data;

import java.io.Serializable;

import etomica.util.Arrays;

/**
 * DataSink that takes a data stream and forwards it across multiple sinks.
 */

/*
 * Sinks are held in an array of DataSinkWrappers (a private inner class defined here).  The
 * wrapper keeps a flag indicating whether a data caster was added between its wrapped sink and this
 * DataFork.  If a caster is added, it is not accessible externally to this class; all public methods deal
 * only with the DataSinks.
 */
public class DataFork implements DataPipeForked, java.io.Serializable {

    /**
     * Constructs with no initial DataSinks.  Sinks can be added and
     * removed after construction using appropriate method calls.  
     */
    public DataFork() {
        tag = new DataTag();
    }
    
    /**
     * Constructs to forward data to the given sinks.  Sinks can be added
     * and removed after construction using appropriate method calls.
     */
    public DataFork(IDataSink[] dataSinks) {
        this();
        setDataSinks(dataSinks);
    }

    public DataTag getTag() {
        return tag;
    }
    
    /**
     * Returns null, indicating that this DataSink can accept any type of Data.
     */
    public DataPipe getDataCaster(IDataInfo incomingDataInfo) {
        return null;
    }
    
    /**
     * Puts the given Data through into all DataSinks. Does
     * nothing if given Data is null.
     */
    public void putData(IData data) {
        if(data != null) {
            for(int i=dataSinkList.length-1; i>=0; i--) {
                dataSinkList[i].dataSink.putData(data);
            }
        }
    }
    
    /**
     * Puts the given DataInfo through into all DataSinks, inserting
     * a data caster before any sinks needing one.
     */
    public void putDataInfo(IDataInfo incomingDataInfo) {
        dataInfo = incomingDataInfo.getFactory().makeDataInfo();
        dataInfo.addTag(tag);
        for(int i=dataSinkList.length-1; i>=0; i--) {
            insertTransformerIfNeeded(i);
            dataSinkList[i].dataSink.putDataInfo(dataInfo);
        }
    }

    public IDataSink[] getDataSinks() {
        IDataSink[] sinks = new IDataSink[dataSinkList.length];
        for (int i=0; i<sinks.length; i++) {
            sinks[i] = dataSinkList[i].trueDataSink;
        }
        return sinks;
    }
    
    /**
     * Sets the given DataSink as the only one connected to this DataFork.  
     * Previously set DataSinks are discarded.  If argument is null, this
     * DataFork will be left with no DataSinks.
     * <p> 
     * Implementation of DataPipe interface.
     * 
     * @param dataSink the new, sole DataSink connected to this DataFork
     */
    public void setDataSink(IDataSink dataSink) {
        setDataSinks(new IDataSink[] {dataSink});
    }

    /* (non-Javadoc)
     * @see etomica.data.DataPipeForked#setDataSinks(etomica.data.DataSink[])
     */
    public void setDataSinks(IDataSink[] dataSinks) {
        dataSinkList = new DataSinkWrapper[0];
        if(dataSinks == null) {
            return;
        }
        for(int i=0; i<dataSinks.length; i++) {
            addDataSink(dataSinks[i]);
        }
    }

    /* (non-Javadoc)
     * @see etomica.data.DataPipeForked#addDataSink(etomica.data.DataSink)
     */
    public void addDataSink(IDataSink newDataSink) {
        if(newDataSink == null) return;
        dataSinkList = (DataSinkWrapper[])Arrays.addObject(dataSinkList, new DataSinkWrapper(newDataSink));
        int numSinks = dataSinkList.length;
        insertTransformerIfNeeded(numSinks-1);
        
        if(dataInfo != null) {
            dataSinkList[numSinks-1].dataSink.putDataInfo(dataInfo);
        }
    }

    /* (non-Javadoc)
     * @see etomica.data.DataPipeForked#removeDataSink(etomica.data.DataSink)
     */
    public void removeDataSink(IDataSink dataSink) {
        for(int i=0; i<dataSinkList.length; i++) {
            if(dataSink == dataSinkList[i].trueDataSink) {
                dataSinkList = (DataSinkWrapper[])Arrays.removeObject(dataSinkList, dataSinkList[i]);
            }
        }
    }
 
    private void insertTransformerIfNeeded(int i) {
        if(dataSinkList[i] == null || dataInfo == null) return;
        //remove transformer if one was previously inserted
        dataSinkList[i].dataSink = dataSinkList[i].trueDataSink;

        DataPipe caster = dataSinkList[i].trueDataSink.getDataCaster(dataInfo);
        if(caster != null) {
            caster.setDataSink(dataSinkList[i].trueDataSink);
            dataSinkList[i].dataSink = caster;
        }
    }

    private static final long serialVersionUID = 1L;
    protected DataSinkWrapper[] dataSinkList = new DataSinkWrapper[0];
    protected IDataInfo dataInfo;
    protected final DataTag tag;
    
    private static class DataSinkWrapper implements Serializable {
        private static final long serialVersionUID = 1L;
        IDataSink dataSink;
        final IDataSink trueDataSink;
        DataSinkWrapper(IDataSink dataSink) {
            this.dataSink = dataSink;
            trueDataSink = dataSink;
        }
    }
}
