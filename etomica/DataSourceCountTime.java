package etomica;

import etomica.Integrator.IntervalEvent;
import etomica.units.Dimension;
import etomica.units.Picosecond;
import etomica.units.Unit;

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

	/**
	 * Identifies the integrators whose steps will be counted. Information
	 * regarding any previously set integrators is discarded.
	 */
	public synchronized void setIntegrator(IntegratorMD[] integrator) {
		//remove existing counters (if any) as listeners to their integrators
		for (int i = 0; i < timer.length; i++) {
			timer[i].integrator.removeIntervalListener(timer[i]);
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
		MyTimer[] newTimer = new MyTimer[timer.length + 1];
		System.arraycopy(timer, 0, newTimer, 0, timer.length);
		newTimer[timer.length] = new MyTimer(integrator);
		timer = newTimer;
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
		MyTimer[] newTimer = new MyTimer[timer.length - 1];
		System.arraycopy(timer, 0, newTimer, 0, i);
		System.arraycopy(timer, i + 1, newTimer, i, newTimer.length - i);
		timer = newTimer;
	}

	private double[] value;

	private MyTimer[] timer = new MyTimer[0];

	//inner class used to handle the counting for each integrator.
	private static class MyTimer implements Integrator.IntervalListener {

		MyTimer(IntegratorMD integrator) {
			this.integrator = integrator;
			integrator.addIntervalListener(this);
		}

		/**
		 * Causes incrementing of timer by current value of evt.getInterval, if
		 * the given event is of type IntervalEvent.INTERVAL (meaning it is not
		 * an event indicating start, stop, etc. of the integrator). If event is
		 * type START, timer is set to zero.
		 */
		public void intervalAction(IntervalEvent evt) {
			if (evt.type() == IntervalEvent.INTERVAL) {
				elapsedTime += evt.getInterval() * integrator.getTimeStep();
			} else if (evt.type() == IntervalEvent.START) {
				elapsedTime = 0.0;
			}
		}

		IntegratorMD integrator;

		double elapsedTime;
	}//end of MyTimer

}