package etomica;

/**
 * Organizer of simulation actions to be executed in series.
 */
public class ActivityGroupSeries extends Activity {
    
    /**
	 * Adds the given action to the list of actions.
	 * If action is already in list of actions to be performed, method returns
	 * without doing anything.
	 */
	public synchronized void addAction(Action newAction) {
		//FIXME this doesn't actually check that newAction isn't already in the array
		actions = (Action[])addObject(actions, newAction);
		numActions++;
	}

    /**
     * Removes the given action from the list of actions performed by this controller.
     * If action is currently running, or is not in list of actions to be performed,
     * method returns without doing anything.
     * @return true if the action was removed
     */
	public synchronized boolean removeAction(Action action) {
		actions = (Action[]) removeObject(actions, action);
		int newNumActions = actions.length;
		if (newNumActions == numActions)
			return false;
		numActions = newNumActions;
		return true;
	}
	
    /**
     * @return a list of the actions yet to be performed by this controller. 
     */
    public synchronized Action[] pendingActions() {return actions;}
    
    /**
     * @return an array containing the action currently being performed  
     * (if there is one).
     */
    public synchronized Action[] currentActions() {return new Action[] {currentAction};}

    /**
     * @return a list of the actions completed by this controller.
     */
    public synchronized Action[] completedActions() {return completedActions;}
    
    /**
     * Causes uncompleted actions added to this controller to be run in sequence.  Should not be
     * executed directly, but instead as part of the Runnable interface it is executed
     * by a thread made upon invoking the start method.
     */
    public void run() {
    	while(numActions > 0) {
    		synchronized(this) {
    			currentAction = actions[0];
    			removeAction(currentAction);
    		}
    		boolean doPause = false;
    		try {
    			currentAction.actionPerformed();
    		}
    		catch (Exception e) {
    			//TODO write message to error stream
    			e.printStackTrace();
    			doPause = true;
    		}
    		doPause = doPause || pauseAfterEachAction;
    		//TODO mark this as whether completed normally
    		synchronized(this) {
    			addAction(currentAction);
    			currentAction = null;
    		}
    		if(doPause || pauseRequested) doWait();
    		if(haltRequested) break;
    	}
    	synchronized(this) {
    		notifyAll();
    	}
    }
            
    /**
     * Method to put controller in a condition of being paused.
     */
    private synchronized void doWait() {
        notifyAll(); //release any threads waiting for pause to take effect
        try {
            wait(); //put in paused state
        } catch(InterruptedException e) {}
        pauseRequested = false;
    }
    
    /**
     * Requests a pause in the performance of the actions. If the current action is
     * an Activity, it is paused; if a simple Action, pause takes effect once it
     * has completed. In either case, calling thread is put in a wait state until 
     * the pause takes effect.
     */
    public synchronized void pause() {
    	if (pauseRequested) return;// already paused
        if(!isActive() || currentAction == null) return;//nothing going on, no need to pause
        pauseRequested = true;
        if(currentAction instanceof Activity) {
        	((Activity)currentAction).pause();//activity enforces pause and has calling thread waits till in effect
        } else {//currentAction is not a pausable activity; put pause in controller loop
	        try {
	            wait();  //make thread requesting pause wait until pause is in effect
	        } catch(InterruptedException e) {}
        }
    }
    
    /**
     * Removes controller from the paused state, resuming execution where it left off.
     */
    public synchronized void unPause() {
    	if (!pauseRequested) return;// not paused
        if(currentAction != null && currentAction instanceof Activity) {
        	((Activity)currentAction).unPause();
        } else {
    		notifyAll();
    	}
    	pauseRequested = false;
    }
    
    /**
     * Queries whether the integrator is in a state of being paused.  This may
     * occur independent of whether the integrator is running or not.  If paused
     * but not running, then pause will take effect upon start.
     */
    public boolean isPaused() {
    	return pauseRequested;
    }
     
    /**
     * @return flag specifying whether controller should pause upon completing each
     * action.
     */
	public boolean isPauseAfterEachActivity() {
		return pauseAfterEachAction;
	}
	/**
	 * @param pauseAfterEachAction specifies whether controller should pause upon
	 * completing each action (true), or if next action should begin immediately
	 * upon completion of current action.
	 */
	public void setPauseAfterEachActivity(boolean pauseAfterEachAction) {
		this.pauseAfterEachAction = pauseAfterEachAction;
	}

    /**
     * Request that the integrator terminate its thread on the next integration step.
     * Does not cause calling thread to wait until this is completed, so it would
     * be prudent to have the calling thread join() to suspend it until the halt
     * is in effect.
     */
    public synchronized void halt() {
        if(!isActive()) return;
        haltRequested = true;
        if(currentAction instanceof Activity) ((Activity)currentAction).halt();
        if(pauseRequested) unPause();
        try {
            wait();  //make thread requesting pause wait until halt is in effect
        } catch(InterruptedException e) {}
    }
    
    protected Action currentAction;
    protected boolean pauseAfterEachAction;
    protected boolean haltRequested;
    protected int numActions;
    protected Action[] actions;
    protected Action[] completedActions;

}//end of Controller


