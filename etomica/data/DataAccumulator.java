/*
 * History
 * Created on Jul 28, 2004 by kofke
 */
package etomica.data;

import etomica.Data;
import etomica.DataSink;
import etomica.DataSource;

/**
 * DataProcessor that accumulates and transforms given data and intermittently transmits
 * processed data to other DataSink(s).
 */
public abstract class DataAccumulator extends DataProcessorForked implements DataSource {

    /**
     * Constructs accumulator with no initial DataSink.
     */
    public DataAccumulator() {
        super();
        setPushInterval(1);
    }
    
	/**
	 * Constructs DataAccumulator with the given DataSink.
	 */
	public DataAccumulator(DataSink dataSink) {
        this();
		setDataSink(dataSink);
	}
    
	/* (non-Javadoc)
	 * @see etomica.Integrator.IntervalListener#intervalAction(etomica.Integrator.IntervalEvent)
	 */
	public Data processData(Data inputData) {
        if(!active) return null;
        addData(inputData);
		if (--putCount == 0) {
		    putCount = pushInterval;
		    return getData();
        }
        return null;
	}

    protected abstract void addData(Data data);
    
    public abstract Data getData();
    
    public abstract void reset();
    
	/**
	 * @return Indicates whether accumulation is being performed in response to interval events.
	 */
	public boolean isActive() {
		return active;
	}
	/**
	 * @param active Sets whether accumulation is to be performed in response to interval events.
	 * Default is true.
	 */
	public final void setActive(boolean active) {
		this.active = active;
	}

    /**
     * Accumulated data are pushed to the data sinks after every pushInterval calls
     * to putData.  This method returns the current value of pushInterval.  Default
     * value is 1.
     */
    public final int getPushInterval() {
        return pushInterval;
    }
    
    /**
     * Accumulated data are pushed to the data sinks after every pushInterval calls
     * to putData.  This method sets the pushInterval, and argument must be greater
     * than zero.  Default value is 1, meaning that every call to addData causes
     * accumulator data to be pushed to its sink.
     * 
     * @throws IllegalArgumentException if argument is less than or equal to zero.
     */
    public final void setPushInterval(int i) {
        if(i > 0) {
            pushInterval = i;
            putCount = pushInterval;
        }
        else throw new IllegalArgumentException("Illegal value of push interval");
    }
    
    /**
     * Overrides the Object class toString method to have it return the output of getName
     * 
     * @return The name given to the phase
     */
    public String toString() {
        return dataInfo.getLabel();
    }

    
    /**
	 * Counter that keeps track of the number of interval events received since last call to updateSums
	 */
	protected int putCount;

    /**
     * Accumulated data are pushed to the data sinks after every pushInterval calls
     */
    private int pushInterval;

	/**
	 * Flag specifying whether the manager responds to integrator events
	 */
	protected boolean active = true;
}
