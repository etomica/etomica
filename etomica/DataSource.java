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
    
    /**
     * Interface for a data source that has associated "x" values.
     */
    public interface X extends DataSource {
        
        public double[] xValues();
        public String getXLabel();
        public etomica.units.Dimension getXDimension();
    }
    
    /**
     * Indicates an object that uses a DataSource.  Useful mainly to the Mediator.
     */
    public interface User {
        public void setDataSource(DataSource source);
        public DataSource getDataSource();
    }
    
    /**
     * Indicates an object that uses multiple DataSources.
     */
    public interface MultiUser {
        public void setDataSources(DataSource[] source);
        public void setDataSources(int i, DataSource source);
        public DataSource[] getDataSources();
        public DataSource getDataSources(int i);
        public void addDataSources(DataSource source);
        public void addDataSources(DataSource[] source);
    }
    
    /**
     * Interface for a class that can provide one or more data source objects, 
     * without necessarily acting as one itself.
     * Example is a Meter, which can provide History and HistogramSimple data sources.
     */
    public interface Wrapper {
        /**
         * An array of strings describing the available data sources from this object.
         */
        public String[] getSourcesAsText();
        
        /**
         * Returns the data source indicated by the argument, which is keyed to one
         * of the strings returned by the getSourcesAsText method.
         */
        public DataSource getDataSource(String text);
    }//end of Wrapper
}//end of DataSource