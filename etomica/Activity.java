/*
 * History
 * Created on Oct 25, 2004 by kofke
 */
package etomica;

/**
 * An Action that supports the capability of pausing/unpausing and terminating
 * on request.
 */
public abstract class Activity implements Action {

	/**
	 * Create class with a simple default label.
	 */
	public Activity() {
		setLabel("Activity");
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see java.lang.Runnable#run()
	 */
	public abstract void run();

	/**
	 * Sets integrator to begin running on its own thread. This is the normal
	 * way to begin the integrator's activity. Fires an event to listeners
	 * indicating that integrator has started, calls the initialize method, and
	 * starts a new thread that then enters the integrators run() method. If
	 * integrator is already running, method call return immediately and has no
	 * effect.
	 */
	public void actionPerformed() {
		running = true;
		haltRequested = false;
		run();
		running = false;
		if (haltRequested) {
			notifyHalt();
			throw new etomica.exception.AbnormalCompletionException();
		}
	}

	private synchronized void notifyHalt() {
		notifyAll();
	}

	protected boolean doContinue() {
		while (pauseRequested)
			doWait();//keep this before resetRequest, since need for reset
		// might naturally follow completion of pause
		//       if(resetRequested) {doReset(); resetRequested = false;}
		if (haltRequested)
			return false;
		//        this.doStep();
		//        if(doSleep) {
		//            try { Thread.sleep(sleepPeriod); }
		//            catch (InterruptedException e) { }
		//        }
		return true;
	}

	/**
	 * Method to put activity in a condition of being paused.
	 */
	private synchronized void doWait() {
		notifyAll(); //release any threads waiting for pause to take effect
		try {
			wait(); //put in paused state
		} catch (InterruptedException e) {
		}
	}

	/**
	 * Requests that the Activity pause its execution. The actual suspension
	 * of execution occurs only after activity notices the pause request.
	 * The calling thread is put in a wait state until the pause takes
	 * effect.
	 */
	public synchronized void pause() {
		if (running && !pauseRequested/* !isPaused */) {
			pauseRequested = true;
			try {
				wait(); //make thread requesting pause wait until pause is in
				// effect
			} catch (InterruptedException e) {
			}
		}
	}

	/**
	 * Removes the integrator from the paused state, resuming execution where it
	 * left off.
	 */
	public synchronized void unPause() {
		pauseRequested = false;
		notifyAll();
	}

	/**
	 * Queries whether the integrator is in a state of being paused. This may
	 * occur independent of whether the integrator is running or not. If paused
	 * but not running, then pause will take effect upon start.
	 */
	public boolean isPaused() {
		return pauseRequested;
	}//isPaused;}

	/**
	 * Indicates if the integrator has been started and has not yet completed.
	 * If so, returns true, even if integrator is presently paused (but not
	 * halted).
	 */
	public boolean isActive() {
		return running;
	}

	/**
	 * Request that the integrator terminate its thread on the next integration
	 * step. Causes calling thread to wait until this is completed.
	 */
	public synchronized void halt() {
		if (!running) return;
		haltRequested = true;
		if (pauseRequested)
			unPause();
		try {
			wait(); //make thread requesting pause wait until halt is in effect
		} catch (InterruptedException e) {
		}
	}

	public String getLabel() {
		return label;
	}

	public void setLabel(String label) {
		this.label = label;
	}

	private boolean running = false;
	private boolean haltRequested = false;
	protected boolean pauseRequested = false;
	private String label;

}