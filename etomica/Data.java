package etomica;


/**
 * TODO To change the template for this generated type comment go to
 * Window - Preferences - Java - Code Style - Code Templates
 *
 * @author David Kofke
 *
 */

/*
 * History
 * Created on Jun 15, 2005 by kofke
 */
public abstract class Data {

    public Data(DataInfo dataInfo) {
        this.dataInfo = dataInfo;
    }
    
    /**
     * Copy constructor, used by subclasses.
     */
    protected Data(Data data) {
        this.dataInfo = new DataInfo(data.dataInfo);
    }
    
    /**
     * @return Returns the dataInfo.
     */
    public DataInfo getDataInfo() {
        return dataInfo;
    }
    
    public abstract Data makeCopy();
    
    public abstract void E(Data data);
    
    protected final DataInfo dataInfo;
}
