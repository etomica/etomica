/*
 * History
 * Created on Jul 28, 2004 by kofke
 */
package etomica;

import java.util.Iterator;
import java.util.LinkedList;

import etomica.Integrator.IntervalEvent;
import etomica.log.DataLogger;

/**
 * Keeps a DataSource and one or more Accumulators and adds to Accumulators
 * in response to IntervalEvents using data from the DataSource.
 */
public class DataManager implements Integrator.IntervalListener {

	private final DataSource dataSource;
	private final LinkedList dataSinkList;
	private Iterator iterator;
	protected DataLogger finishDataDumper;
    private int priority;
	
	/**
	 * Constructs DataManager with the given DataSource and
	 * Accumulators.
	 */
	public DataManager(DataSource dataSource, DataSink[] dataSinks) {
		if(dataSource == null) throw new NullPointerException("Error: cannot construct accumulator manager without a data source");
		this.dataSource = dataSource;
		dataSinkList = new LinkedList(); 
		setDataSinks(dataSinks);
		setUpdateInterval(1);
		setEventType(Integrator.IntervalEvent.INTERVAL);
        setPriority(200);
	}
	
	/**
	 * Constructor with AccumulatorAverage as the default Accumulator.
	 * @param dataSource
	 */
	public DataManager(DataSource dataSource) {
		this(dataSource, new DataSink[] {new AccumulatorAverage()});
	}

	public void setDumpOnFinishFile(String fileName) {
		finishDataDumper = new DataLogger(fileName);
	}
	
	/* (non-Javadoc)
	 * @see etomica.Integrator.IntervalListener#intervalAction(etomica.Integrator.IntervalEvent)
	 */
	public void intervalAction(IntervalEvent evt) {
	    if(!active) return;
	    if(evt.type() == eventType || (evt.type() == Integrator.IntervalEvent.DONE && finishDataDumper != null)) {
	    	if (eventType == Integrator.IntervalEvent.INTERVAL) {
	    		if (--iieCount != 0) return;
                iieCount = updateInterval;
	    	}
            if (evt.type() == eventType) {
                double[] data = dataSource.getData();
                iterator = dataSinkList.iterator();
                while(iterator.hasNext()) {
                    ((DataSink)iterator.next()).add(data);
                }
            }
		    if(evt.type() == Integrator.IntervalEvent.DONE && finishDataDumper != null) {
                iterator = dataSinkList.iterator();
		        while(iterator.hasNext()) {
		        	DataSink dataSink = (DataSink)iterator.next();
		        	if (dataSink instanceof DataSource) {
		        		if (dataSink instanceof AccumulatorAverage) {
		        			finishDataDumper.add(((AccumulatorAverage)dataSink).getData(AccumulatorAverage.AVERAGE));
		        			finishDataDumper.add(((AccumulatorAverage)dataSink).getData(AccumulatorAverage.STANDARD_DEVIATION));
                            finishDataDumper.add(((AccumulatorAverage)dataSink).getData(AccumulatorAverage.ERROR));
		        		}
		        		else {
		        			finishDataDumper.add(((DataSource)dataSink).getData());
		        		}
		        	}
		        }
		    }
	    }
	}
	
	/**
	 * @return Indicates if accumulation is being performed in response to interval events.
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
	 * @return Returns the accumulators.
	 */
	public DataSink[] getDataSinks() {
		return (DataSink[])dataSinkList.toArray();
	}
	/**
	 * @param dataSinks The accumulators to set.
	 */
	public void setDataSinks(DataSink[] dataSinks) {
		dataSinkList.clear();
		for(int i=0; i<dataSinks.length; i++) {
			dataSinkList.add(dataSinks[i]);
			if (dataSinks[i] instanceof Accumulator)
				((Accumulator)dataSinks[i]).setDataSourceDimension(dataSource.getDimension());
		}
	}
	
	public void addDataSink(DataSink dataSink) {
		dataSinkList.add(dataSink);
		if (dataSink instanceof Accumulator)
			((Accumulator)dataSink).setDataSourceDimension(dataSource.getDimension());
	}
	
	/**
	 * Removes the specified accumulator from this manager.
     * @param accumulator accumulator to be removed from this list, if present.
     * @return <tt>true</tt> if the manager contained the specified accumulator.
     */
	public boolean removeDataSink(DataSink dataSink) {
		return dataSinkList.remove(dataSink);
	}
	
    /**
     * Accessor method for the updateInterval
     * @see #updateInterval
     */
    public final int getUpdateInterval() {return updateInterval;}
    /**
     * Accessor method for the updateInterval.
     * Sets to given value and resets count of interval events
     * @see #updateInterval
     */
    public final void setUpdateInterval(int i) {
        if(i > 0) {
            updateInterval = i;
            iieCount = updateInterval;
        }
        else throw new IllegalArgumentException("Illegal value of update interval");
    }
	
    /**
     * Invokes the reset() method of all accumulators held by this manager.
     */
    public void resetAccumulators() { 
        iterator = dataSinkList.iterator();
        while(iterator.hasNext()) {
            DataSink dataSink = (DataSink)iterator.next();
            if (dataSink instanceof Accumulator) {
                ((Accumulator)dataSink).reset();
            }
        } 
    } 
    
    public final void setEventType(Integrator.IntervalEvent.Type type) {
    	eventType = type;
    }
    
    /**
     * determines the type of IntervalEvent the data manager is interested in
     */
    protected Integrator.IntervalEvent.Type eventType;

    /**
     * @return Returns the interval-listener priority.
     */
    public int getPriority() {
        return priority;
    }
    /**
     * Sets the interval-listener priority.  Default value is 200, which
     * puts this after central-image enforcement.
     * @param priority The priority to set.
     */
    public void setPriority(int priority) {
        this.priority = priority;
    }
    /**
	 * Counter that keeps track of the number of interval events received since last call to updateSums
	 */
	protected int iieCount;
	/**
	 * A string describing the property measured by the meter
	 */
	protected String label = "Property";
	/**
	 * Number of integration interval events received before another call to updateSums
	 */
	protected int updateInterval;

	/**
	 * Flag specifying whether the meter responds to integrator events
	 * If false, the meter does not perform regular measurements or keep averages
	 * In this situation the meter is probably measuring a property for use by some other object
	 * Default is <code>true</code> for a Meter, but <code>false</code> for a MeterFunction.
	 */
	protected boolean active=true;
}
