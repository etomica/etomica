package etomica.data;

import etomica.Integrator;
import etomica.integrator.IntegratorHard;
import etomica.units.Count;
import etomica.units.Dimension;
import etomica.units.Unit;

/**
 * This is a data source to count the number of collisions processed by a
 * hard-potential integrator.
 */

public class DataSourceCountCollisions extends DataSourceAdapter {

	/**
	 * Sets up data source to count collisions for the given integrator.
	 */
	public DataSourceCountCollisions(IntegratorHard integrator) {
		this();
		addIntegrator(integrator);
	}

	/**
	 * Sets up data source with no integrator specified.  Requires
	 * call to addIntegrator or setIntegrator before use.
	 */
	public DataSourceCountCollisions() {
		super(Dimension.QUANTITY);
		setLabel("Number of Collisions");
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

	public Unit defaultIOUnit() {
		return Count.UNIT;
	}

	/**
	 * Returns the number of steps performed by each of the integrators tracked
	 * by this class. Each value corresponds to the integrators given in
	 * setIntegrator and/or addIntegrator (with most recently added integrator
	 * given last in returned array)
	 */
	public double[] getData() {
		for (int i = 0; i < counter.length; i++) {
			value[i] = (double) counter[i].count;
		}
		return value;
	}
    
    public int getDataLength() {
        return counter.length;
    }

	/**
	 * Identifies the integrators whose steps will be counted. Information
	 * regarding any previously set integrators is discarded.
	 * 
	 * @param integrator
	 */
	public synchronized void setIntegrator(IntegratorHard[] integrator) {
		//remove existing counters (if any) as listeners to their integrators
		for (int i = 0; i < counter.length; i++) {
			counter[i].integrator.removeCollisionListener(counter[i]);
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
	public synchronized void addIntegrator(IntegratorHard integrator) {
		//check that integrator isn't already added
		for (int i = 0; i < counter.length; i++) {
			if (counter[i].integrator == integrator)
				return;
		}
		MyCounter[] newCounter = new MyCounter[counter.length + 1];
		System.arraycopy(counter, 0, newCounter, 0, counter.length);
		newCounter[counter.length] = new MyCounter(integrator);
		counter = newCounter;
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
		MyCounter[] newCounter = new MyCounter[counter.length - 1];
		System.arraycopy(counter, 0, newCounter, 0, i);
		System.arraycopy(counter, i + 1, newCounter, i, newCounter.length - i);
		counter = newCounter;
	}

	private double[] value;

	private MyCounter[] counter = new MyCounter[0];

	//inner class used to handle the counting for each integrator.
	private static class MyCounter implements IntegratorHard.CollisionListener {

		MyCounter(IntegratorHard integrator) {
			this.integrator = integrator;
			integrator.addCollisionListener(this);
		}

		/**
		 * Implements CollisionListener. Adds one to the counter with each
		 * collision.
		 */
		public void collisionAction(IntegratorHard.Agent agent) {
			count++;
		}

		IntegratorHard integrator;

		int count;
	}//end of MyCounter

}

