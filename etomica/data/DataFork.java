package etomica.data;

import java.io.Serializable;

import etomica.Data;
import etomica.DataInfo;
import etomica.DataSink;
import etomica.utility.Arrays;

/**
 * DataSink that takes a data stream and forwards it across multiple sinks.
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
        

    public DataProcessor getDataCaster(DataInfo dataInfo) {
        return null;
    }
    
    /**
     * Method called by subclasses to move data into sinks.
     */
    public void putData(Data data) {
        if(data != null) {
            for(int i=dataSinkList.length-1; i>=0; i--) {
                dataSinkList[i].dataSink.putData(data);
            }
        }
    }
    
    public void putDataInfo(DataInfo dataInfo) {
        this.dataInfo = dataInfo;
        for(int i=dataSinkList.length-1; i>=0; i--) {
            insertTransformerIfNeeded(i);
            dataSinkList[i].dataSink.putDataInfo(dataInfo);
        }
    }

    /**
     * @return the current set of data sinks.
     */
    public DataSink getDataSink(int i) {
        return dataSinkList[i].dataSink;
    }
    
    /**
     * Sets the given DataSink as the only one connected to this DataFork.  
     * Previously set DataSinks are discarded.  If argument is null, this
     * DataFork will be left with no DataSinks.<p> 
     * Implementation of DataPipe interface.
     * 
     * @param dataSink the new, sole DataSink connected to this DataFork
     */
    public void setDataSink(DataSink dataSink) {
        setDataSinks(new DataSink[] {dataSink});
    }

    /**
     * Sets the list of DataSinks that receive the Data entering this DataFork.
     * All previously added DataSinks are discarded.
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
        insertTransformerIfNeeded(dataSinkList.length-1);
        
        if(dataInfo != null) {
            newDataSink.putDataInfo(dataInfo);
        }
    }

    /**
     * Removes the specified data sink.
     * 
     * @param dataSink data sink to be removed from this list, if present.
     */
    public void removeDataSink(DataSink dataSink) {
        for(int i=0; i<dataSinkList.length; i++) {
            if(dataSinkList[i].dataSink == dataSink) {
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
