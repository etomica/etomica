package etomica.data;

import java.io.Serializable;
import java.util.HashMap;

import etomica.data.types.DataGroup;
import etomica.util.Arrays;

/**
 * Receives data from multiple streams and organizes it in the form of a table,
 * with each stream of data forming a column of the table. Fires event when data
 * in table is changed; can be configured to fire when any column of data
 * changes, or to perform action only when all columns are changed.
 * 
 * @author David Kofke
 *  
 */

/*
 * History Created on Apr 9, 2005 by kofke
 */
public class DataSet implements DataSink, Serializable {

    public DataSet() {
        this(null);
    }
    
    public DataSet(DataCasterJudge dataCasterJudge) {
        allData = new Data[0];
        this.dataCasterJudge = dataCasterJudge;
    }
    
    public DataCasterJudge getDataCasterJudge() {
        return dataCasterJudge;
    }
    
    /* (non-Javadoc)
     * @see etomica.DataSink#getDataCaster(etomica.DataInfo)
     */
    public DataProcessor getDataCaster(DataInfo dataInfo) {
        if (dataCasterJudge == null) {
            return null;
        }
        return dataCasterJudge.getDataCaster(dataInfo);
    }

    /* (non-Javadoc)
     * @see etomica.DataSink#putData(etomica.Data)
     */
    public void putData(Data data) {
        boolean newData = false;
        Integer indexObj = (Integer)dataIndexHash.get(data);
        if (indexObj == null) {
            newData = true;
            addNewData(data);
            indexObj = new Integer(casterIndex++);
            if ((casterIndex+31)/32 > casterChangedBits.length) {
                casterChangedBits = Arrays.resizeArray(casterChangedBits,casterChangedBits.length+1);
            }
            casterChangedBits[casterChangedBits.length-1] |= 1<<((casterIndex-1)%32);
            dataIndexHash.put(data,indexObj);
        }
        int index = indexObj.intValue();
        // unset the bit to indicate this data object has been updated
        casterChangedBits[index >> 5] &= (0xFFFFFFFF ^ (1<<(index & 0xFFFFFFFF)));
        if (updatingOnAnyChange) {
            fireDataChangedEvent();
        }
        else {
            int l = casterChangedBits.length;
            // check to see if all data has been updated
            // always redraw when new data comes in
            if (!newData) {
                for (int i=0; i<l; i++) {
                    if (casterChangedBits[i] != 0) {
                        return;
                    }
                }
            }
            fireDataChangedEvent();
            // reset all bits but the last one
            for (int i=0; i<l-1; i++) {
                casterChangedBits[i] = 0xFFFFFFFF;
            }
            if (casterIndex%32 == 0) {
                // reset the last one
                casterChangedBits[l-1] = 0xFFFFFFFF;
            }
            else {
                // partially reset the last one
                casterChangedBits[l-1] = (1<<(casterIndex%32)) - 1;
            }
        }
    }

    protected void addNewData(Data newData) {
        if (newData instanceof DataGroup) {
            int oldSize = allData.length;
            allData = (Data[])Arrays.resizeArray(allData, oldSize+((DataGroup)newData).getNData());
            for (int i=oldSize; i<allData.length; i++) {
                allData[i] = ((DataGroup)newData).getData(i-oldSize);
            }
        }
        else {
            allData = (Data[])Arrays.addObject(allData, newData);
        }
        fireDataCountChangedEvent();
    }
    
    /* (non-Javadoc)
     * @see etomica.DataSink#putDataInfo(etomica.DataInfo)
     */
    public void putDataInfo(DataInfo dataInfo) {
    }
    
    public void reset() {
        dataIndexHash.clear();
        allData = new Data[0];
        casterIndex = 0;
        casterChangedBits = new int[0];
        fireDataCountChangedEvent();
    }
    
    
    public Data getData(int i) {
        return allData[i];
    }

    public int getDataCount() {
        return allData.length;
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

    protected DataSetListener[] listeners = new DataSetListener[0];
    private boolean updatingOnAnyChange;
    private int casterIndex;
    private int[] casterChangedBits = new int[0];
    private HashMap dataIndexHash = new HashMap();
    private Data[] allData;
    private final DataCasterJudge dataCasterJudge;
    
    public interface DataCasterJudge {
        public DataProcessor getDataCaster(DataInfo inputDataInfo);
    }
}
