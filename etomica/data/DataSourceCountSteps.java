package etomica.data;

import etomica.EtomicaElement;
import etomica.EtomicaInfo;
import etomica.Integrator;
import etomica.Integrator.IntervalEvent;
import etomica.units.Count;
import etomica.units.Dimension;
import etomica.units.Unit;
import etomica.utility.Arrays;

/**
 * Data source that keeps track of the number of steps performed by an
 * integrator. More precisely, sum the integrator's interval value each time the
 * integrator fires an INTERVAL event. Normally, this will equal the number of
 * times the integrator's doStep method has been called. A START event from the
 * integrator will reset the count.
 */

public final class DataSourceCountSteps extends DataSourceAdapter implements
		EtomicaElement {

	/**
	 * Sets up data source with no integrator specified. Requires call to
	 * addIntegrator or setIntegrator before use.
	 */
	public DataSourceCountSteps() {
		super(Dimension.QUANTITY);
		setLabel("Integrator steps");
	}

	/**
	 * Sets up data source to count steps for the given integrator.
	 */
	public DataSourceCountSteps(Integrator integrator) {
		this();
		addIntegrator(integrator);
	}

	public static EtomicaInfo getEtomicaInfo() {
		EtomicaInfo info = new EtomicaInfo(
				"Records the number of steps performed by the integrator");
		return info;
	}

	/**
	 * @return Count.UNIT
	 */
	public Unit defaultIOUnit() {
		return Count.UNIT;
	}

	/**
	 * Resets all counters to zero
	 */
	public void reset() {
		for (int i = 0; i < counter.length; i++) {
			counter[i].count = 0;
		}
	}

	/**
	 * Resets to zero the counter for the given integrator.
	 */
	public void reset(Integrator integrator) {
		for (int i = 0; i < counter.length; i++) {
			if (counter[i].integrator == integrator) {
				counter[i].count = 0;
				return;
			}
		}
	}

	/**
	 * Returns the number of steps performed by each of the integrators tracked
	 * by this class. Each value corresponds to the integrators given in
	 * setIntegrator and/or addIntegrator (with most recently added integrator
	 * given last in returned array)
	 */
	public double[] getData() {
		for (int i = 0; i < counter.length; i++) {
			value[i] = counter[i].count;
		}
		return value;
	}

	/**
	 * Identifies the integrators whose steps will be counted. Information
	 * regarding any previously set integrators is discarded.
	 * 
	 * @param integrator
	 */
	public synchronized void setIntegrator(Integrator[] integrator) {
		//remove existing counters (if any) as listeners to their integrators
		for (int i = 0; i < counter.length; i++) {
			counter[i].integrator.removeIntervalListener(counter[i]);
		}
		//make new counters for the integrators
		counter = new MyCounter[integrator.length];
		for (int i = 0; i < counter.length; i++) {
			counter[i] = new MyCounter(integrator[i]);
		}
		value = new double[integrator.length];
	}

	/**
	 * Adds the given integrator to those having steps counted.
	 * 
	 * @param integrator
	 */
	public synchronized void addIntegrator(Integrator integrator) {
		//check that integrator isn't already added
		for (int i = 0; i < counter.length; i++) {
			if (counter[i].integrator == integrator)
				return;
		}
        counter = (MyCounter[])Arrays.addObject(counter, new MyCounter(integrator));
        value = new double[counter.length];
	}

	/**
	 * Removes the given integrator from those having steps counted. No action
	 * is performed if given integrator was not previously added.
	 * 
	 * @param integrator
	 *            to be removed.
	 */
	public synchronized void removeIntegrator(Integrator integrator) {
		//check that integrator isn't already added
		int i;
		for (i = 0; i < counter.length; i++) {
			if (counter[i].integrator == integrator)
				break;
		}
		if (i == counter.length)
			return; //didn't find it
        counter = (MyCounter[])Arrays.removeObject(counter, counter[i]);
        value = new double[counter.length];
	}

	private double[] value;

	private MyCounter[] counter = new MyCounter[0];

	//inner class used to handle the counting for each integrator.
	private static class MyCounter implements Integrator.IntervalListener {

		MyCounter(Integrator integrator) {
			this.integrator = integrator;
			integrator.addIntervalListener(this);
		}
        
        /**
         * Priority is 1, which ensures that counters are updated before
         * any meters might be called to use them.
         */
        public int getPriority() {return 1;}

		/**
		 * Causes incrementing of counter by current value of evt.getInterval,
		 * if the given event is of type IntervalEvent.INTERVAL (meaning it is
		 * not an event indicating start, stop, etc. of the integrator). If
		 * event is type START, counter is set to zero.
		 */
		public void intervalAction(IntervalEvent evt) {
			if (evt.type() == IntervalEvent.INTERVAL) {
				count += evt.getInterval();
			} else if (evt.type() == IntervalEvent.START) {
				count = 0;
			}
		}

		Integrator integrator;

		int count;
	}//end of MyCounter

}