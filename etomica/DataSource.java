package etomica;

/**
 * Interface for an object that can provide data 
 * (in the form of an array of double) for plotting or other analysis.
 */
 
public interface DataSource {

	public double[] getData();
    
    /**
     * Returns a label used to describe the data when presented, 
     * for example as a plot legend.
     */
    public String getLabel();
    
    /**
     * Returns the physical dimensions (e.g., length) of the data.
     */
    public etomica.units.Dimension getDimension();
    
    /**
     * Returns the DataTranslator that converts this source's data
     * into an appropriate object. 
     */
    public DataTranslator getTranslator();
        
}//end of DataSource