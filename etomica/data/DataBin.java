package etomica.data;

import etomica.Data;
import etomica.DataInfo;
import etomica.DataSink;
import etomica.DataSource;

/**
 * DataSink that receives and holds data without pushing it downstream. May
 * notify a DataBinManager when the data is changed, which can then access the
 * data via the DataSource methods implemented by this class. Dimension (e.g.,
 * length, time) is set at construction and cannot be changed (this is because
 * the class does not push its data, but instead serves as a data source to its
 * DataBinManager; it therefore does not have a way to notify downstream sinks
 * if its dimension field changes).
 * 
 * @author David Kofke
 *  
 */

/*
 * History Created on Apr 9, 2005 by kofke
 */
public class DataBin implements DataSink, DataSource, java.io.Serializable {

    /**
     * Constructs DataBin with null DataBinManager
     */
    public DataBin(DataInfo info) {
        this(info, null);
    }

    /**
     * Construct new DataBin that will notify the given DataBinManager (via the
     * manager's dataChangeNotify method) (if not null) any time the bin's
     * putData method is invoked.
     * 
     * @param dataBinManager
     *            manger of this bin; declared final; null value is permitted
     */
    public DataBin(DataInfo info, DataBinManager dataBinManager) {
        this.dataInfo = info;
        this.dataBinManager = dataBinManager;
    }

    /**
     * Stores the given array of values. Data are held in a copy array, so any
     * subsequent changes to given array will not affect values held by the bin.
     */
    public void putData(Data newData) {
        if(data == null) {
            data = newData.makeCopy();
        } else {
            data.E(newData);
        }
        dataInfo = data.getDataInfo();
        if (dataBinManager != null) {
            dataChanged = true;
            dataBinManager.dataChangeNotify(this);
        }
    }

    /**
     * Returns the array holding the data most recently given to putData.
     * Returns a zero-length array if putData was not previously invoked.
     */
    public Data getData() {
        return data;
    }
    
    public DataInfo getDataInfo() {
        return dataInfo;
    }

    private Data data;
    private DataInfo dataInfo;
    private final DataBinManager dataBinManager;

    //flag used by DataBinManager to keep track of bin changes
    boolean dataChanged = true;
}
