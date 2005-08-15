package etomica.data;

import java.io.Serializable;

import etomica.utility.Arrays;

/**
 * DataSink that takes a data stream and forwards it across multiple sinks.
 */

/*
 * Sinks are held in an array of DataSinkWrappers (a private inner class defined here).  The
 * wrapper keeps a flag indicating whether a data caster was added between its wrapped sink and this
 * DataFork.  If a caster is added, it is not accessible externally to this class; all public methods deal
 * only with the DataSinks.
 */

/*
 * Created July 21, 2005
 */
public class DataFork implements DataPipe, java.io.Serializable {

    /**
     * Constructs with no initial DataSinks.  Sinks can be added and
     * removed after construction using appropriate method calls.  
     */
    public DataFork() {
    }
    
    /**
     * Constructs to forward data to the given sinks.  Sinks can be added
     * and removed after construction using appropriate method calls.
     */
    public DataFork(DataSink[] dataSinks) {
        this();
        setDataSinks(dataSinks);
    }
        

    /**
     * Returns null, indicating that this DataSink can accept any type of Data.
     */
    public DataProcessor getDataCaster(DataInfo dataInfo) {
        return null;
    }
    
    /**
     * Puts the given Data through into all DataSinks. Does
     * nothing if given Data is null.
     */
    public void putData(Data data) {
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
    public void putDataInfo(DataInfo dataInfo) {
        this.dataInfo = dataInfo;
        for(int i=dataSinkList.length-1; i>=0; i--) {
            insertTransformerIfNeeded(i);
            dataSinkList[i].dataSink.putDataInfo(dataInfo);
        }
    }

    /**
     * @return the i-th DataSink
     */
    public DataSink getDataSink(int i) {
        return dataSinkList[i].dataSink;
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
    public void setDataSink(DataSink dataSink) {
        setDataSinks(new DataSink[] {dataSink});
    }

    /**
     * Sets the list of DataSinks that receive the Data entering this DataFork.
     * All previously added DataSinks are discarded.  If argument is null, all
     * existing DataSinks are discarded and none are added.
     * 
     * @param dataSinks The data sinks to set.
     */
    public void setDataSinks(DataSink[] dataSinks) {
        dataSinkList = new DataSinkWrapper[0];
        if(dataSinks == null) {
            return;
        }
        for(int i=0; i<dataSinks.length; i++) {
            addDataSink(dataSinks[i]);
        }
    }

    /**
     * Adds the given DataSink to those receiving the Data entering this DataFork,
     * keeping all previously entered DataSinks.  If argument is null, no action
     * is performed.
     * 
     * @param dataSink
     */
    public void addDataSink(DataSink newDataSink) {
        if(newDataSink == null) return;
        dataSinkList = (DataSinkWrapper[])Arrays.addObject(dataSinkList, new DataSinkWrapper(newDataSink));
        int numSinks = dataSinkList.length;
        insertTransformerIfNeeded(numSinks-1);
        
        if(dataInfo != null) {
            dataSinkList[numSinks-1].dataSink.putDataInfo(dataInfo);
        }
    }

    /**
     * Removes the specified data sink.  Does nothing if the given DataSink is
     * not currently a DataSink for this DataFork.
     * 
     * @param dataSink data sink to be removed from this list, if present.
     */
    public void removeDataSink(DataSink dataSink) {
        for(int i=0; i<dataSinkList.length; i++) {
            DataSink testSink = null;
            //if we inserted a transformer, we have to look past it for the sink
            if(dataSinkList[i].insertedTransformer) {
                DataProcessor dataCaster = (DataProcessor)dataSinkList[i].dataSink;
                testSink = dataCaster.dataSink;
            //otherwise, just check the dataSink
            } else {
                testSink = dataSinkList[i].dataSink;
            }
            if(testSink == dataSink) {
                dataSinkList = (DataSinkWrapper[])Arrays.removeObject(dataSinkList, dataSinkList[i]);
            }
        }
    }
 
    private void insertTransformerIfNeeded(int i) {
        if(dataSinkList[i] == null || dataInfo == null) return;
        //remove transformer if one was previously inserted
        if(dataSinkList[i].insertedTransformer) {
            DataProcessor dataCaster = (DataProcessor)dataSinkList[i].dataSink;
            dataSinkList[i] = new DataSinkWrapper(dataCaster.dataSink);
        }
        DataProcessor caster = dataSinkList[i].dataSink.getDataCaster(dataInfo);
        if(caster != null) {
            caster.setDataSink(dataSinkList[i].dataSink);
            dataSinkList[i] = new DataSinkWrapper(caster);
            dataSinkList[i].insertedTransformer = true;
        }
    }

    protected DataSinkWrapper[] dataSinkList = new DataSinkWrapper[0];
    protected DataInfo dataInfo;
    
    private static class DataSinkWrapper implements Serializable {
        final DataSink dataSink;
        boolean insertedTransformer = false;
        DataSinkWrapper(DataSink dataSink) {
            this.dataSink = dataSink;
        }
    }

}
