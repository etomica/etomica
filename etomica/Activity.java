/*
 * History
 * Created on Oct 25, 2004 by kofke
 */
package etomica;

/**
 * @author kofke
 *
 * TODO To change the template for this generated type comment go to
 * Window - Preferences - Java - Code Style - Code Templates
 * 
 * TODO update comments to reflect this is not an integrator
 */
public abstract class Activity implements Action {

	/**
	 * 
	 */
	public Activity() {
		setLabel("Activity");
	}

	/* (non-Javadoc)
	 * @see java.lang.Runnable#run()
	 */
	public abstract void run();
	
   /**
     * Sets integrator to begin running on its own thread.  This is the normal
     * way to begin the integrator's activity.  Fires an event to listeners indicating
     * that integrator has started, calls the initialize method, and starts a new
     * thread that then enters the integrators run() method.  If integrator is already
     * running, method call return immediately and has no effect.
     */
    public void actionPerformed() {
        running = true;
        haltRequested = false;
        run();
        if (haltRequested) {
        	notifyHalt();
        	throw new etomica.exception.AbnormalCompletionException();
        }
    }
    
    private synchronized void notifyHalt() {
    	notifyAll();
    }
    
    protected boolean doContinue() {
        while(pauseRequested) doWait();//keep this before resetRequest, since need for reset might naturally follow completion of pause
 //       if(resetRequested) {doReset(); resetRequested = false;}
        if(haltRequested) return false;
//        this.doStep();
//        if(doSleep) {
//            try { Thread.sleep(sleepPeriod); }
//            catch (InterruptedException e) { }
//        }
        return true;
    }
    
    /**
     * Method to put integrator in a condition of being paused.
     */
    private synchronized void doWait() {
//        isPaused = true;
		//System.out.println("pausing");
        notifyAll(); //release any threads waiting for pause to take effect
        try {
            wait(); //put in paused state
        } catch(InterruptedException e) {}
//        isPaused = false;
		//System.out.println("done pausing");
    }
    
    //suspend and resume functions
    /**
     * Requests that the integrator pause its execution.  The actual suspension
     * of execution occurs only after completion of the current integration step.
     * The calling thread is put in a wait state until the pause takes effect.
     */
    public synchronized void pause() {
        if(running && !pauseRequested/*!isPaused*/) {
            pauseRequested = true;
            try {
                wait();  //make thread requesting pause wait until pause is in effect
            } catch(InterruptedException e) {}
        }
    }
    /**
     * Removes the integrator from the paused state, resuming execution where it left off.
     */
    public synchronized void unPause() {pauseRequested = false; notifyAll();}
    /**
     * Queries whether the integrator is in a state of being paused.  This may
     * occur independent of whether the integrator is running or not.  If paused
     * but not running, then pause will take effect upon start.
     */
    public boolean isPaused() {return pauseRequested;}//isPaused;}
    
    /**
     * Indicates if the integrator has been started and has not yet completed.
     * If so, returns true, even if integrator is presently paused (but not halted).
     */
    public boolean isActive() {return running;}
    
    //stop function
    //consider having calling thread here join() or wait() for halt to take effect
    /**
     * Request that the integrator terminate its thread on the next integration step.
     * Does not cause calling thread to wait until this is completed, so it would
     * be prudent to have the calling thread join() to suspend it until the halt
     * is in effect.
     */
    public synchronized void halt() {
        if(running) haltRequested = true;
        if(pauseRequested) unPause();
        try {
            wait();  //make thread requesting pause wait until halt is in effect
        } catch(InterruptedException e) {}
    }

    public String getLabel() {
    	return label;
    }
    public void setLabel(String label) {
    	this.label = label;
    }
    
	transient public Thread runner;
	private boolean running = false;
	private boolean haltRequested = false;
	private boolean resetRequested = false;
	 //looking to eliminate isPaused field and used just pauseRequested
	//	  private boolean isPaused = true;
	protected boolean pauseRequested = false;
	private String label;

}
