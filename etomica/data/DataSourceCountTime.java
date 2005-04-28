package etomica.data;

import etomica.EtomicaElement;
import etomica.EtomicaInfo;
import etomica.Integrator;
import etomica.IntegratorNonintervalEvent;
import etomica.IntegratorIntervalEvent;
import etomica.IntegratorIntervalListener;
import etomica.IntegratorNonintervalListener;
import etomica.integrator.IntegratorMD;
import etomica.units.Dimension;
import etomica.units.Picosecond;
import etomica.units.Unit;
import etomica.utility.Arrays;

/**
 * Data source that keeps track of the elapsed simulation time of an MD
 * integrator. More precisely, sums the integrator's interval value times the
 * time-step each time the integrator fires an INTERVAL event. A START event
 * from the integrator will reset the elapsedTime.
 */

public final class DataSourceCountTime extends DataSourceAdapter implements
		EtomicaElement {

	/**
	 * Sets up data source with no integrator specified.  Requires
	 * call to addIntegrator or setIntegrator before use.
	 */
	public DataSourceCountTime() {
		super(Dimension.TIME);
		setLabel("Simulation Time");
	}
	/**
	 * Sets up data source to count time for the given integrator.
	 */
	public DataSourceCountTime(IntegratorMD integrator) {
		this();
		addIntegrator(integrator);
	}

	public static EtomicaInfo getEtomicaInfo() {
		EtomicaInfo info = new EtomicaInfo(
				"Records the elapsed simulation time for an MD simulation");
		return info;
	}

	/**
	 * @return Picosecond.UNIT
	 */
	public Unit defaultIOUnit() {
		return Picosecond.UNIT;
	}

	/**
	 * Resets all timers to zero
	 */
	public void reset() {
		for (int i = 0; i < timer.length; i++) {
			timer[i].elapsedTime = 0.0;
		}
	}

	/**
	 * Resets to zero the timer for the given integrator.
	 */
	public void reset(Integrator integrator) {
		for (int i = 0; i < timer.length; i++) {
			if (timer[i].integrator == integrator) {
				timer[i].elapsedTime = 0.0;
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
		for (int i = 0; i < timer.length; i++) {
			value[i] = timer[i].elapsedTime;
		}
		return value;
	}
    
    public int getDataLength() {
        return timer.length;
    }

	/**
	 * Identifies the integrators whose steps will be counted. Information
	 * regarding any previously set integrators is discarded.
	 */
	public synchronized void setIntegrator(IntegratorMD[] integrator) {
		//remove existing counters (if any) as listeners to their integrators
		for (int i = 0; i < timer.length; i++) {
			timer[i].integrator.removeListener(timer[i]);
		}
		//make new counters for the integrators
		timer = new MyTimer[integrator.length];
		for (int i = 0; i < timer.length; i++) {
			timer[i] = new MyTimer(integrator[i]);
		}
		value = new double[integrator.length];
	}

	/**
	 * Adds the given integrator to those having steps counted.
	 */
	public synchronized void addIntegrator(IntegratorMD integrator) {
		//check that integrator isn't already added
		for (int i = 0; i < timer.length; i++) {
			if (timer[i].integrator == integrator)
				return;
		}
        timer = (MyTimer[])Arrays.addObject(timer, new MyTimer(integrator));
        value = new double[timer.length];
	}

	/**
	 * Removes the given integrator from those having steps counted. No action
	 * is performed if given integrator was not previously added.
	 * 
	 * @param integrator
	 *            integrator to be removed.
	 */
	public synchronized void removeIntegrator(Integrator integrator) {
		//check that integrator isn't already added
		int i;
		for (i = 0; i < timer.length; i++) {
			if (timer[i].integrator == integrator)
				break;
		}
		if (i == timer.length)
			return; //didn't find it
        timer = (MyTimer[])Arrays.removeObject(timer, timer[i]);
        value = new double[timer.length];
	}

	private double[] value;

	private MyTimer[] timer = new MyTimer[0];

	//inner class used to handle the counting for each integrator.
	private static class MyTimer implements IntegratorNonintervalListener, IntegratorIntervalListener {

		MyTimer(IntegratorMD integrator) {
			this.integrator = integrator;
			integrator.addListener(this);
		}

		/**
		 * Causes incrementing of timer by current value of evt.getInterval, if
		 * the given event is of type IntervalEvent.INTERVAL (meaning it is not
		 * an event indicating start, stop, etc. of the integrator). If event is
		 * type START, timer is set to zero.
		 */
		public void intervalAction(IntegratorIntervalEvent evt) {
			elapsedTime += evt.getInterval() * integrator.getTimeStep();
        }
        
        public void nonintervalAction(IntegratorNonintervalEvent evt) {
			if(evt.type() == IntegratorIntervalEvent.START) {
				elapsedTime = 0.0;
			}
		}
        /**
         * Priority is 1, which ensures that counters are updated before
         * any meters might be called to use them.
         */
        public int getPriority() {return 1;}


		IntegratorMD integrator;

		double elapsedTime;
	}//end of MyTimer

}