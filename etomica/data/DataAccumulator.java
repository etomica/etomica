/*
 * History
 * Created on Jul 28, 2004 by kofke
 */
package etomica.data;

import etomica.DataPipe;
import etomica.DataSink;
import etomica.DataSource;

/**
 * DataPipe that accumulates and transforms given data and intermittently transmits
 * processed data to other DataSink(s).
 */
public abstract class DataAccumulator extends DataPipe implements DataSource {

    public DataAccumulator() {
        setPushInterval(1);
    }
    
	/**
	 * Constructs DataAccumulator with the given DataSinks.
	 */
	public DataAccumulator(DataSink[] dataSinks) {
        this();
		setDataSinks(dataSinks);
	}
    
    /**
     * Constructs DataAccumulator with a single DataSink.
     */
    public DataAccumulator(DataSink dataSink) {
        this(new DataSink[] {dataSink});
    }
	
	/* (non-Javadoc)
	 * @see etomica.Integrator.IntervalListener#intervalAction(etomica.Integrator.IntervalEvent)
	 */
	public void putData(double[] newData) {
        if(!active) return;
        addData(newData);
		if (--putCount == 0) {
		    putCount = pushInterval;
		    pushData(getData());
        }
	}

    protected abstract void addData(double[] data);
    
    public abstract double[] getData();
    
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
     * to putData.  This method returns the current value of pushInterval.
     */
    public final int getPushInterval() {
        return pushInterval;
    }
    
    /**
     * Accumulated data are pushed to the data sinks after every pushInterval calls
     * to putData.  This method sets the pushInterval, and argument must be greater
     * than zero.
     */
    public final void setPushInterval(int i) {
        if(i > 0) {
            pushInterval = i;
            putCount = pushInterval;
        }
        else throw new IllegalArgumentException("Illegal value of push interval");
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
	protected boolean active=true;
    
}
