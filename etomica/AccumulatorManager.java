/*
 * History
 * Created on Jul 28, 2004 by kofke
 */
package etomica;

import etomica.Integrator.IntervalEvent;
import etomica.utility.java2.LinkedList;
import etomica.utility.java2.Iterator;

/**
 * @author kofke
 *
 * Keeps a DataSource and one or more Accumulators and adds to Accumulators
 * in response to IntervalEvents using data from the DataSource.
 */
public class AccumulatorManager implements Integrator.IntervalListener {

	private final DataSource dataSource;
	private final LinkedList accumulatorList;
	private final Iterator iterator;
	
	/**
	 * Constructs AccumulatorManager with the given DataSource and
	 * Accumulators.
	 */
	public AccumulatorManager(DataSource dataSource, Accumulator[] accumulators) {
		this.dataSource = dataSource;
		accumulatorList = new LinkedList(); 
		iterator = accumulatorList.iterator();
		setAccumulators(accumulators);
	}
	
	/**
	 * Constructor with AccumulatorAverage as the default Accumulator.
	 * @param dataSource
	 */
	public AccumulatorManager(DataSource dataSource) {
		this(dataSource, new Accumulator[] {new AccumulatorAverage()});
	}

	/* (non-Javadoc)
	 * @see etomica.Integrator.IntervalListener#intervalAction(etomica.Integrator.IntervalEvent)
	 */
	public void intervalAction(IntervalEvent evt) {
	    if(!active) return;
	    if(evt.type() != Integrator.IntervalEvent.INTERVAL) return; //don't act on start, done, initialize events
	    if(--iieCount == 0) {
	        iieCount = updateInterval;
	        double[] data = dataSource.getData();
	        iterator.reset();
	        while(iterator.hasNext()) {
	        	((Accumulator)iterator.next()).add(data);
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
	public Accumulator[] getAccumulators() {
		return (Accumulator[])accumulatorList.toArray();
	}
	/**
	 * @param accumulators The accumulators to set.
	 */
	public void setAccumulators(Accumulator[] accumulators) {
		accumulatorList.clear();
		for(int i=0; i<accumulators.length; i++) {
			accumulatorList.add(accumulators[i]);
		}
	}
	
	public void addAccumulator(Accumulator accumulator) {
		accumulatorList.add(accumulator);
	}
	
	/**
	 * Removes the specified accumulator from this manager.
     * @param accumulator accumulator to be removed from this list, if present.
     * @return <tt>true</tt> if the manager contained the specified accumulator.
     */
	public boolean removeAccumulator(Accumulator accumulator) {
		return accumulatorList.remove(accumulator);
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
	protected boolean active;
}
