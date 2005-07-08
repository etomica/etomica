package etomica.data;

import etomica.Data;
import etomica.DataSink;
import etomica.utility.Arrays;

/**
 * An object that may hold one or more data sinks to which it can 
 * push its data.
 */
public abstract class DataPusher implements java.io.Serializable {

    /**
     * Method called by subclasses to move data into sinks.
     */
    protected void pushData(Data data) {
        for(int i=dataSinkList.length-1; i>=0; i--) {
            dataSinkList[i].putData(data);
        }
    }

    /**
     * @return Returns the data sinks.
     */
    public DataSink[] getDataSinks() {
        return dataSinkList;
    }

    /**
     * @param dataSinks The data sinks to set.
     */
    public void setDataSinks(DataSink[] dataSinks) {
        if(dataSinks == null) {
            dataSinkList = new DataSink[0];
            return;
        }
        dataSinkList = (DataSink[])dataSinks.clone();
    }

    public void addDataSink(DataSink dataSink) {
        if(dataSink == null) return;
        dataSinkList = (DataSink[])Arrays.addObject(dataSinkList, dataSink);
    }

    /**
     * Removes the specified data sink from this manager.
     * @param dataSink data sink to be removed from this list, if present.
     * @return <tt>true</tt> if the manager contained the specified data sink.
     */
    public void removeDataSink(DataSink dataSink) {
        dataSinkList = (DataSink[])Arrays.removeObject(dataSinkList, dataSink);
    }
    
    protected DataSink[] dataSinkList = new DataSink[0];


}
