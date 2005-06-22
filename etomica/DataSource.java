package etomica;

/**
 * Interface for an object that can provide data 
 * (in the form of an array of doubles) for plotting or other analysis.
 */
 
public interface DataSource {

    /**
     * @return the data given by this source
     */
	public Data getData();
    
    public DataInfo getDataInfo();
    
}//end of DataSource