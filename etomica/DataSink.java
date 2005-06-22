package etomica;


/**
 * A consumer of data.  Data goes in and might (or might not) come out.
 */
public interface DataSink {
    
    /**
     * Gives data to DataSink for processing, display, or whatever it does.
     */
	public abstract void putData(Data data);
    
}