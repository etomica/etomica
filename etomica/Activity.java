/*
 * History
 * Created on Oct 25, 2004 by kofke
 */
package etomica;

import etomica.utility.NameMaker;

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
        setName(NameMaker.makeName(this.getClass()));
	}

	/**
	 * Method defining the behavior of the activity.  Implementation should
	 * ensure regular checking of doContinue() to permit any requests to pause
	 * or halt to be put in effect.
	 */
	public abstract void run();

	/**
	 * Sets integrator to begin isActive on its own thread. This is the normal
	 * way to begin the integrator's activity. Fires an event to listeners
	 * indicating that integrator has started, calls the initialize method, and
	 * starts a new thread that then enters the integrators run() method. If
	 * integrator is already isActive, method call return immediately and has no
	 * effect.
	 */
	public void actionPerformed() {
        synchronized(this) {
            haltRequested = false;
            pauseRequested = false;
            isPaused = false;
            isActive = true;
        }
		run();
		isActive = false;
		if (haltRequested) {
			synchronized(this) {notifyAll();}//release thread waiting for halt to take effect
			throw new etomica.exception.AbnormalCompletionException();
		}
	}

	protected synchronized boolean doContinue() {
//        System.out.println(this+" in doContinue");
		while (pauseRequested)
			doWait();//keep this before resetRequest, since need for reset
		// might naturally follow completion of pause
		//       if(resetRequested) {doReset(); resetRequested = false;}
//        System.out.println(this+" got past doWait in doContinue");
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
	 * Should be called in run() method of subclasses when they find
	 * that pauseRequested is true, or if they wish to pause the execution thread
	 * for their own purposes.
	 */
	protected synchronized void doWait() {
		notifyAll(); //release any threads waiting for pause to take effect
//        System.out.println(this+"in doWait, setting isPaused to true");
		isPaused = true;
		try {
			wait(); //put in paused state
		} catch (InterruptedException e) {}
        System.out.println(this+"in doWait, setting isPaused to false");
		isPaused = false;
	}

	/**
	 * Requests that the Activity pause its execution. The actual suspension of
	 * execution occurs only after activity notices the pause request. The
	 * calling thread is put in a wait state until the pause takes effect.
	 */
	public synchronized void pause() {
//        System.out.println("in Activity.pause "+isPaused+" "+isActive());
		if (isPaused || !isActive()) return;// already paused or not active
		pauseRequested = true;
		try {
			// make thread requesting pause wait until pause is in effect
//            System.out.println("in Activity.pause waiting");
			wait();
//            System.out.println("in Activity.pause waited");
		} catch (InterruptedException e) {
			//interrupted before pause took effect; abort pause request
			pauseRequested = false;
		}
	}

    /**
     * Removes activity from the paused state, resuming execution where it left off.
     */
    public synchronized void unPause() {
//        System.out.println(this+" in unPause "+isPaused+" "+isActive());
    	if (!isPaused || !isActive()) return;// not currently paused or not active
    	pauseRequested = false;
    	notifyAll();
    }

	/**
	 * Request that the activity terminate as soon as safely possible. 
	 * Causes calling thread to wait until the halt is in effect.
	 */
    public synchronized void halt() {
        if(!isActive()) return;
        haltRequested = true;
        unPause();//in case currently paused
        try {
            wait();  //make thread requesting pause wait until halt is in effect
        } catch(InterruptedException e) {}
    }

	/**
	 * Queries whether the integrator is in a state of being paused. This may
	 * occur independent of whether the integrator is isActive or not. If paused
	 * but not isActive, then pause will take effect upon start.
	 */
	public synchronized boolean isPaused() {
		return isPaused;
	}

	/**
	 * Indicates if the integrator has been started and has not yet completed.
	 * If so, returns true, even if integrator is presently paused (but not
	 * halted).
	 */
	public synchronized boolean isActive() {
		return isActive;
	}


	public String getLabel() {
		return label;
	}

	public void setLabel(String label) {
		this.label = label;
	}

    /**
     * @return Returns the name.
     */
    public String getName() {
        return name;
    }
    /**
     * @param name The name to set.
     */
    public void setName(String name) {
        this.name = name;
    }
	private boolean isActive = false;

	protected boolean haltRequested = false;

	protected boolean pauseRequested = false;

	protected boolean isPaused = false;

	private String label;
    
    private String name;

}